has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_solver sdpt3;
cvx_precision low;
cvx_clear;
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general parameters
T = 24*7; % planning period length, in hours
IP = 0.25; % idle power consumption of a server
PP = 0.5; % peak power consumption of a server
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
OP = 1.3; % over-provisioning ratio of power and server

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
IM = [20]; % mean of interactive workload consumption, in kW
for i = 1:1:IN
    a(i,:) = trace_process(char(IL(i)), T, 12, 4, 1,IM(i));
end
au = [0.3]; % average utilization of interactive workloads
con = 0.8; % consolidation ratio, i.e., percentage of spare capacity that can be shared by batch jobs

% batch job
BN = 40; % number of batch jobs
BM = 60; % mean of batch job consumption, in kW
[A,BS,S,E] = batch_job_generator(T,BN,'Uniform',3.99,23.99,'Uniform',1,4,BM); % generate batch jobs based on statistical properties

% renewable generation
RN = 1; % number of renewable generation
RL = {'traces\solar-one-week.csv'};
R = zeros(IN,T);
RP = [130]; % peak of renewable generation
for i = 1:1:RN
    R(i,:) = trace_process(char(RL(i)), T, 12, 1, 2, RP(i));
end

% cooling efficiency
PUE = 1.5*ones(1,T);
PUE_low = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 150, 1000];
GP = 0.12*ones(1,T); % grid electricity price
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cost_flat,capacity_flat] = no_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,P,GP,2,'results\flat');
[cost_PV,capacity_PV] = PV_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,R(1,:),P,GP,BP,1,'results\PV');
[cost_supply,capacity_supply] = supply_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,R(1,:),P,GP,RP,BP,1,'results\supply');
[cost_demand,capacity_demand] = demand_integration(T,IP,PP,OP,CO2_grid,a,au,con,BS,A,S,E,PUE_low,R(1,:),P,GP,RP,BP,1,'results\demand',capacity_PV);
[cost_opt,capacity_opt] = opt(T,IP,PP,OP,CO2_grid,a,au,con,BS,A,S,E,PUE_low,R(1,:),P,GP,RP,BP,2,'results\opt',capacity_flat(4));

%[cost,b,EU,status,C,EC,RC] = tco_opt(T,a,BS,A,R,PUE,P_IT,P_ES,P_R,Loss,1,'results\tco1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
bar([cost_flat;cost_PV;cost_supply;cost_demand;cost_opt],0.2,'stacked');
ylabel('annual expenditure ($)');
legend('Power','Cooling','PV', 'IT','Electricity')
set(gca,'xticklabel',{'Grid-only','PV','supply-side','demand-side','our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\expenditure_comp1.eps');
saveas(gcf,'results\expenditure_comp1.fig');

figure;
bar([capacity_flat;capacity_PV;capacity_supply;capacity_demand;capacity_opt]);
ylabel('capacity (kW)');
legend('Power','Cooling','PV', 'IT')
set(gca,'xticklabel',{'Grid-only','PV','supply-side','demand-side','our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\capacity_comp1.eps');
saveas(gcf,'results\capacity_comp1.fig');
%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}