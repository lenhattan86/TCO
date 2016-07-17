has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_solver sedumi;
cvx_precision low;
cvx_clear;
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general parameters
T = 24*7; % planning period length, in hours
IP = 0.1; % idle power consumption of a server, in kW
PP = 0.4; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
DC_power = 500; % Data center power capacity, in kW
OP = 0.6; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.8; % maximum utilization after consolidation

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
au = 0.3; % average utilization of interactive workloads
IPU = 0.8; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 3; % peak to mean ratio, should be smaller than current PMR
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
end

% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.25; % batch job ratio, compared with interactive workload
[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties
bu = 0.9; % average utilization of batch job when running alone.

% uncontrolled renewable generation, e.g., PV, wind
RN = 1; % number of renewable generation
RL = {'traces\solar-one-week-07012012.csv'};
R = zeros(IN,T);
RP = 0.3; % capacity factor
RR = [0, 500]; % renewable capacity range, in kW
for i = 1:1:RN
    R(i,:) = trace_process(char(RL(i)), T, 12, 1, 1, RP(i));
end
RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
RE = 32; % emmision, in gram CO2 eq

% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
CRC = [1000/20+SR*200/4, 0.01+0.07]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
CRE = 443; % emmision, in gram CO2 eq

% cooling efficiency
PUE = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
GP = 0.08*ones(1,T); % grid electricity price
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
discount_factor = 0;
max_payback = 20;
payback_plot = 0;
net_zero = 1;
option = 3; % 1 for grid & PV, 2 for PV & gas, 3 for gas & grid
start1 = 0.02;
end1 = 0.30;
step1 = 0.02;
start2 = 0.10;
end2 = 0.30;
step2 = 0.02;
start3 = 0.01;
end3 = 0.30;
step3 = 0.05;

% default setting
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.07];GP = 0.08*ones(1,T);CRR=[0,1400];
if option == 1 % varying grid power price & PV capacity factor
    for i = 1:1:(floor((end1-start1)/step1)+1)
        % grid only
        grid(i) = start1 + (i-1)*step1;
        [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R,[0,0],RC,RE,SR,[0,0],CRC,CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
        yearly_cost_grid = sum(cost_grid([2,4,5]));
        for j = 1:1:(floor((end2-start2)/step2)+1)  
            % supply only
            PV(j) = start2 + (j-1)*step2;
            [cost_supply,capacity_supply,emission_supply] = houston_grid(0,net_zero,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R/RP*PV(j),RR,RC,RE,SR,CRR,CRC,CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
            yearly_cost_supply = sum(cost_supply([2,4,5]));
            infra_cost_supply = sum(cost_supply([1,3]))*max_payback;
            [payback_supply(i,j)] = payback_cal(infra_cost_supply,(yearly_cost_grid-yearly_cost_supply)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
            % integrated supply and demand opt
            [cost_int,capacity_int,emission_int] = houston_grid(1,net_zero,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R/RP*PV(j),RR,RC,RE,SR,CRR,CRC,CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
            yearly_cost_int = sum(cost_int([2,4,5]));
            infra_cost_int = sum(cost_int([1,3]))*max_payback;
            [payback_int(i,j)] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
        end
    end
end
% default setting
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.07];GP = 0.08*ones(1,T);CRR=[0,1400];
if option == 2 % varying PV capacity factor and gas price
    % grid only
    [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R,[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_grid = sum(cost_grid([2,4,5]));
    for i = 1:1:(floor((end2-start2)/step2)+1)
        PV(i) = start2 + (i-1)*step2;
        for j = 1:1:(floor((end3-start3)/step3)+1) 
            gas(j) = start3 + (j-1)*step3;
            % supply only
            [cost_supply,capacity_supply,emission_supply] = houston_grid(0,net_zero,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R/RP*PV(i),RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas(j)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
            yearly_cost_supply = sum(cost_supply([2,4,5]));
            infra_cost_supply = sum(cost_supply([1,3]))*max_payback;
            [payback_supply(i,j)] = payback_cal(infra_cost_supply,(yearly_cost_grid-yearly_cost_supply)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
            % integrated supply and demand opt
            [cost_int,capacity_int,emission_int] = houston_grid(1,net_zero,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R/RP*PV(i),RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas(j)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
            yearly_cost_int = sum(cost_int([2,4,5]));
            infra_cost_int = sum(cost_int([1,3]))*max_payback;
            [payback_int(i,j)] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
        end
    end
end
% default setting
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.07];GP = 0.08*ones(1,T);CRR=[0,1400];
if option == 3 % varying grid power price and gas price
    for i = 1:1:(floor((end1-start1)/step1)+1)
        % grid only
        grid(i) = start1 + (i-1)*step1;
        [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R,[0,0],RC,RE,SR,[0,0],CRC,CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
        yearly_cost_grid = sum(cost_grid([2,4,5]));
        for j = 1:1:(floor((end3-start3)/step3)+1)  
            % supply only
            gas(j) = start3 + (j-1)*step3;
            [cost_supply,capacity_supply,emission_supply] = houston_grid(0,net_zero,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R,RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas(j)],CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
            yearly_cost_supply = sum(cost_supply([2,4,5]));
            infra_cost_supply = sum(cost_supply([1,3]))*max_payback;
            [payback_supply(i,j)] = payback_cal(infra_cost_supply,(yearly_cost_grid-yearly_cost_supply)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
            % integrated supply and demand opt
            [cost_int,capacity_int,emission_int] = houston_grid(1,net_zero,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R,RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas(j)],CRE,P,grid(i)*ones(1,T),RP,BP,plot_option,'results\tmp');
            yearly_cost_int = sum(cost_int([2,4,5]));
            infra_cost_int = sum(cost_int([1,3]))*max_payback;
            [payback_int(i,j)] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if option == 1
    surf(grid,PV,payback_supply');
    hold on;
    surf(grid,PV,payback_int');
    colormap(cool);
    xlim([start1,end1]);
    ylim([start2,end2]);
    xlabel('grid power price ($/kWh)')
    ylabel('PV capacity factor')
end
if option == 2
    surf(PV,gas,payback_supply');
    hold on;
    surf(PV,gas,payback_int');
    colormap(cool);
    xlim([start2,end2]);
    ylim([start3,end3]);
    xlabel('PV capacity factor')
    ylabel('gas price ($/kWh)')
end
if option == 3
    surf(grid,gas,payback_supply');
    hold on;
    surf(grid,gas,payback_int');
    colormap(cool);
    xlim([start1,end1]);
    ylim([start3,end3]);
    xlabel('grid power price ($/kWh)')
    ylabel('gas price ($/kWh)')
end
save('results\tco_payback_surface.mat');

%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}