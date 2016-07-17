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
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
DC_power = 500; % Data center power capacity, in kW
OP = 0.8; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.9; % maximum utilization after consolidation

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
au = 0.4; % average utilization of interactive workloads
IPU = 0.85; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 3; % peak to mean ratio, should be smaller than current PMR
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
end

% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.25; % batch job ratio, compared with interactive workload
[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties
bu = 1; % average utilization of batch job when running alone.

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
CRC = [1000/20+SR*200/4, 0.005+0.06]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
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


%IT
ITR = [0,DC_power*OP];
ITC = 4000/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 1;
discount_factor = 0;
max_payback = 20;
payback_plot = 0;
% grid only
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.06];GP = 0.08*ones(1,T);CRR=[0,1400];
[cost_grid,capacity_grid,emission_grid] = houston_grid_IT(0,0,0,0,T,IP,PP,DC_power,[DC_power*OP,DC_power*OP],ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_grid = sum(cost_grid([2,4,5]));

% net-zero, supply only
[cost_nz,capacity_nz,emission_nz] = houston_grid_IT(0,1,0,0,T,IP,PP,DC_power,ITR,ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_nz = sum(cost_nz([2,4,5]));
infra_cost_nz = sum(cost_nz([1,3]))*max_payback;
[payback_nz,npv_nz] = payback_cal(infra_cost_nz,(yearly_cost_grid-yearly_cost_nz)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, demand only
[cost_ds,capacity_ds,emission_ds] = houston_grid_IT(1,1,1,1,T,IP,PP,DC_power,[DC_power*OP,DC_power*OP],ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_ds = sum(cost_ds([2,4,5]));
infra_cost_ds = sum(cost_ds([1,3]))*max_payback;
[payback_ds,npv_ds] = payback_cal(infra_cost_ds,(yearly_cost_grid-yearly_cost_ds)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, integrated supply and demand opt
[cost_int,capacity_int,emission_int] = houston_grid_IT(1,1,1,1,T,IP,PP,DC_power,ITR,ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_int = sum(cost_int([2,4,5]));
infra_cost_int = sum(cost_int([1,3]))*max_payback;
[payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
bar([cost_grid;cost_nz;cost_ds;cost_int],0.2,'stacked')
ylabel('annual expenditure ($)');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill','IT capital cost','Location','Northeast');
set(gca,'xticklabel',{'Grid','Supply','Demand','Our solution'});

figure;
bar([capacity_grid;capacity_nz;capacity_ds;capacity_int])
legend('PV','GE','Grid','IT');
ylabel('capacity (kW)');
set(gca,'xticklabel',{'Grid','Supply','Demand','Our solution'});

%{
figure;
bar([cost_nz(1),cost_nz(3),0,cost_nz(6);cost_nz(2),cost_nz(4),cost_nz(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_ds(1),cost_ds(3),0,cost_ds(6);cost_ds(2),cost_ds(4),cost_ds(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_int(1),cost_int(3),0,cost_int(6);cost_int(2),cost_int(4),cost_int(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');
%}
save('results\houston_payback_IT.mat');