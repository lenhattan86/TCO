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
OP = 1; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.8; % maximum utilization after consolidation

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
au = 0.3; % average utilization of interactive workloads
IPU = 0.9; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 3; % peak to mean ratio, should be smaller than current PMR
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
end

% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.5; % batch job ratio, compared with interactive workload
[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',11.99,23.99,'Uniform',1,2,BM/(1-BM)*sum(mean(a,2)./au)); % generate batch jobs based on statistical properties

% uncontrolled renewable generation, e.g., PV, wind
RN = 1; % number of renewable generation
RL = {'traces\solar-one-week-1.csv'};
R = zeros(IN,T);
RP = 0.2; % capacity factor
RR = [0, 500]; % renewable capacity range, in kW
for i = 1:1:RN
    R(i,:) = trace_process(char(RL(i)), T, 12, 1, 1, RP(i));
end
RC = [3200/20, 0.01+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]

% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
CRC = [1000/20+SR*200/4, 0.02+0.08]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]

% cooling efficiency
PUE = 1.3*ones(1,T);
PUE_low = 1.3*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
GP = 0.12*ones(1,T); % grid electricity price
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
gas = 0.02+[0.04, 0.08, 0.12];
for i = 1:1:3
    [cost_houston(i,:),capacity_houston(i,:),G] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),RR,RC,SR,CRR,[CRC(1),gas(i)],P,GP,RP,BP,plot_option,'results\tmp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
bar([cost_houston],0.2,'stacked');
ylabel('annual expenditure ($)');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Location','NorthWest');
%ylim([0,sum(cost_houston_givenBoth)*1.3]);
set(gca,'xticklabel',{'0.04','0.08','0.12'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\expenditure_comp_houston_price_bar.eps');
saveas(gcf,'results\expenditure_comp_houston_price_bar.fig');

figure;
bar([capacity_houston(1,:);capacity_houston(2,:);capacity_houston(3,:)]);
ylabel('capacity (kW)');
legend('PV','GE','Location','NorthWest');
%ylim([0,max(capacity_houston_givenBoth)*1.3]);
set(gca,'xticklabel',{'0.04','0.08','0.12'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\capacity_comp_houston_price_bar.eps');
saveas(gcf,'results\capacity_comp_houston_price_bar.fig');

save('results\houston_price_bar.mat');
