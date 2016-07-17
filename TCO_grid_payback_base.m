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
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general parameters
T = 24*7; % planning period length, in hours
IP = 0.1; % idle power consumption of a server, in kW
PP = 0.4; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
DC_power = 300; % Data center power capacity, in kW
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
plot_option = 1;
discount_factor = 0;
max_payback = 20;
payback_plot = 1;
% grid only
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.07];GP = 0.08*ones(1,T);CRR=[0,1400];
[cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_grid = sum(cost_grid([2,4,5]));
% net-zero, supply only
[cost_nz,capacity_nz,emission_nz] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_nz = sum(cost_nz([2,4,5]));
infra_cost_nz = sum(cost_nz([1,3]))*max_payback;
[payback_nz,npv_nz] = payback_cal(infra_cost_nz,(yearly_cost_grid-yearly_cost_nz)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% net-zero, base supply option
[cost_real,capacity_real,emission_real] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_real = sum(cost_real([2,4,5]));
infra_cost_real = sum(cost_real([1,3]))*max_payback;
[payback_real,npv_real] = payback_cal(infra_cost_real,(yearly_cost_grid-yearly_cost_real)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% net-zero, base supply option, demand opt
[cost_realopt,capacity_realopt,emission_realopt] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_realopt = sum(cost_realopt([2,4,5]));
infra_cost_realopt = sum(cost_realopt([1,3]))*max_payback;
[payback_realopt,npv_realopt] = payback_cal(infra_cost_realopt,(yearly_cost_grid-yearly_cost_realopt)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% net-zero, integrated supply and demand opt
[cost_int,capacity_int,emission_int] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_int = sum(cost_int([2,4,5]));
infra_cost_int = sum(cost_int([1,3]))*max_payback;
[payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% net-zero, PV only
[cost_PV,capacity_PV,emission_PV] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_PV = sum(cost_PV([2,4,5]));
infra_cost_PV = sum(cost_PV([1,3]))*max_payback;
[payback_PV,npv_PV] = payback_cal(infra_cost_PV,(yearly_cost_grid-yearly_cost_PV)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% net-zero, PV only, demand opt
[cost_PVopt,capacity_PVopt,emission_PVopt] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_PVopt = sum(cost_PVopt([2,4,5]));
infra_cost_PVopt = sum(cost_PVopt([1,3]))*max_payback;
[payback_PVopt,npv_PVopt] = payback_cal(infra_cost_PVopt,(yearly_cost_grid-yearly_cost_PVopt)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% not net-zero, supply only
[cost_notsupply,capacity_notsupply,emission_notsupply] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_notsupply = sum(cost_notsupply([2,4,5]));
infra_cost_notsupply = sum(cost_notsupply([1,3]))*max_payback;
[payback_notsupply,npv_notsupply] = payback_cal(infra_cost_notsupply,(yearly_cost_grid-yearly_cost_notsupply)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
% not net-zero, integrated supply and demand opt
[cost_notint,capacity_notint,emission_notint] = houston_grid(1,0,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
yearly_cost_notint = sum(cost_notint([2,4,5]));
infra_cost_notint = sum(cost_notint([1,3]))*max_payback;
[payback_notint,npv_notint] = payback_cal(infra_cost_notint,(yearly_cost_grid-yearly_cost_notint)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
bar([cost_grid;cost_nz;cost_real;cost_int],0.2,'stacked');
ylabel('annual expenditure ($)');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill');
ylim([0,sum(capacity_real)*1.3]);
set(gca,'xticklabel',{'Grid only','Supply opt','Supply given','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\expenditure_comp_pb.eps');
saveas(gcf,'results\expenditure_comp_pb.fig');

figure;
bar([capacity_grid;capacity_nz;capacity_real;capacity_int]);
ylabel('capacity (kW)');
legend('PV','GE');
ylim([0,max(capacity_real)*1.3]);
set(gca,'xticklabel',{'Grid only','Supply opt','Supply given','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\capacity_comp_pb.eps');
saveas(gcf,'results\capacity_comp_pb.fig');

figure;
bar([cost_nz(1),cost_nz(3),0;cost_nz(2),cost_nz(4),cost_nz(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_real(1),cost_real(3),0;cost_real(2),cost_real(4),cost_real(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_realopt(1),cost_realopt(3),0;cost_realopt(2),cost_realopt(4),cost_realopt(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_int(1),cost_int(3),0;cost_int(2),cost_int(4),cost_int(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_PV(1),cost_PV(3),0;cost_PV(2),cost_PV(4),cost_PV(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_PVopt(1),cost_PVopt(3),0;cost_PVopt(2),cost_PVopt(4),cost_PVopt(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_notsupply(1),cost_notsupply(3),0;cost_notsupply(2),cost_notsupply(4),cost_notsupply(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_notint(1),cost_notint(3),0;cost_notint(2),cost_notint(4),cost_notint(5)]*max_payback)
legend('PV','GE','Grid');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

save('results\houston_payback.mat');

%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}