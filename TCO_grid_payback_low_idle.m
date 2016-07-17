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
IP = 0.04; % idle power consumption of a server, in kW
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 2;
discount_factor = 0;
max_payback = 20;
payback_plot = 1;

% grid only
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.06];GP = 0.08*ones(1,T);CRR=[0,1400];
[cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_grid = sum(cost_grid([2,4,5]));
total_cost_grid = yearly_cost_grid*max_payback;

% net-zero, supply only
[cost_nz,capacity_nz,emission_nz] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_nz = sum(cost_nz([2,4,5]));
infra_cost_nz = sum(cost_nz([1,3]))*max_payback;

figure;
bar([cost_nz(1),cost_nz(3),0;cost_nz(2),cost_nz(4),cost_nz(5)]*max_payback)
legend('PV','GE','Grid','Location','Northwest');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_nz = yearly_cost_nz*max_payback + infra_cost_nz;
cost_reduction_nz = 1- total_cost_nz/total_cost_grid

[payback_nz,npv_nz] = payback_cal(infra_cost_nz,(yearly_cost_grid-yearly_cost_nz)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, integrated supply and demand opt
[cost_int,capacity_int,emission_int] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_int = sum(cost_int([2,4,5]));
infra_cost_int = sum(cost_int([1,3]))*max_payback;

figure;
bar([cost_int(1),cost_int(3),0;cost_int(2),cost_int(4),cost_int(5)]*max_payback)
legend('PV','GE','Grid');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_int = yearly_cost_int*max_payback + infra_cost_int;
cost_reduction_int = 1- total_cost_int/total_cost_grid

[payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)


% net-zero, supply only, low gas price
[cost_supplycheapgas,capacity_supplycheapgas,emission_supplycheapgas] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[237,237],RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+0.03],CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_supplycheapgas = sum(cost_supplycheapgas([2,4,5]));
infra_cost_supplycheapgas = sum(cost_supplycheapgas([1,3]))*max_payback;

figure;
bar([cost_supplycheapgas(1),cost_supplycheapgas(3),0;cost_supplycheapgas(2),cost_supplycheapgas(4),cost_supplycheapgas(5)]*max_payback)
legend('PV','GE','Grid');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_supplycheapgas = yearly_cost_supplycheapgas*max_payback + infra_cost_supplycheapgas;
cost_reduction_supplycheapgas = 1- total_cost_supplycheapgas/total_cost_grid

[payback_supplycheapgas,npv_supplycheapgas] = payback_cal(infra_cost_supplycheapgas,(yearly_cost_grid-yearly_cost_supplycheapgas)*ones(1,max_payback),discount_factor,max_payback,payback_plot)


% net-zero, supply & demand, low gas price
[cost_intcheapgas,capacity_intcheapgas,emission_intcheapgas] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[249,249],RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+0.03],CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_intcheapgas
yearly_cost_intcheapgas = sum(cost_intcheapgas([2,4,5]));
infra_cost_intcheapgas = sum(cost_intcheapgas([1,3]))*max_payback;

figure;
bar([cost_intcheapgas(1),cost_intcheapgas(3),0;cost_intcheapgas(2),cost_intcheapgas(4),cost_intcheapgas(5)]*max_payback)
legend('PV','GE','Grid');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_intcheapgas = yearly_cost_intcheapgas*max_payback + infra_cost_intcheapgas;
cost_reduction_intcheapgas = 1- total_cost_intcheapgas/total_cost_grid

[payback_intcheapgas,npv_intcheapgas] = payback_cal(infra_cost_intcheapgas,(yearly_cost_grid-yearly_cost_intcheapgas)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, supply given
[cost_supplygiven,capacity_supplygiven,emission_supplygiven] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_supplygiven
yearly_cost_supplygiven = sum(cost_supplygiven([2,4,5]));
infra_cost_supplygiven = sum(cost_supplygiven([1,3]))*max_payback;

figure;
bar([cost_supplygiven(1),cost_supplygiven(3),0;cost_supplygiven(2),cost_supplygiven(4),cost_supplygiven(5)]*max_payback)
legend('PV','GE','Grid','Location','Northwest');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_supplygiven = yearly_cost_supplygiven*max_payback + infra_cost_supplygiven;
cost_reduction_supplygiven = 1- total_cost_supplygiven/total_cost_grid

[payback_supplygiven,npv_supplygiven] = payback_cal(infra_cost_supplygiven,(yearly_cost_grid-yearly_cost_supplygiven)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, supply given & demand
[cost_intgiven,capacity_intgiven,emission_intgiven] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_intgiven
yearly_cost_intgiven = sum(cost_intgiven([2,4,5]));
infra_cost_intgiven = sum(cost_intgiven([1,3]))*max_payback;

figure;
bar([cost_intgiven(1),cost_intgiven(3),0;cost_intgiven(2),cost_intgiven(4),cost_intgiven(5)]*max_payback)
legend('PV','GE','Grid');
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

total_cost_intgiven = yearly_cost_intgiven*max_payback + infra_cost_intgiven;
cost_reduction_intgiven = 1- total_cost_intgiven/total_cost_grid

[payback_intgiven,npv_intgiven] = payback_cal(infra_cost_intgiven,(yearly_cost_grid-yearly_cost_intgiven)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save('results\houston_payback.mat');
