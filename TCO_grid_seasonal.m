close all; clear all;
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
TS = 24; % number of slots per day
ND = 12; % number of days
T = TS*ND; % planning period length, in hours
IP = 0.04; % idle power consumption of a server, in kW
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
DC_power = 500; % Data center power capacity, in kW
OP = 0.8; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.9; % maximum utilization after consolidation
plot_PV = 0; % plot PV traces
plot_demand = 0; % plot power demand

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces/SAPnew/sapTrace2.tab'};
a = zeros(IN,T);
au = 0.4; % average utilization of interactive workloads
IPU = 0.7; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 2.5; % peak to mean ratio, should be smaller than current PMR
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
RL = {'traces/solar-one-week-07012012.csv'};
RCF = load('traces/PV_CF.csv');
R = zeros(IN,T);
RP = 0.3; % capacity factor
RR = [0, 500]; % renewable capacity range, in kW
for i = 1:1:RN
    R_raw = trace_process(char(RL(i)), TS*7, 12, 1, 4, TS);
    R(i,:) = reshape(R_raw'/mean(R_raw) * RCF,1,T)/100;
end
RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
RE = 32; % emmision, in gram CO2 eq

% plot PV traces
if plot_PV
    figure;
    plot(1/TS:1/TS:ND, R(1,:)*150,'r-','LineWidth',2)
    ylabel('generation (kW)')
    ylim([0,150])
    set(gca,'XTick',[0.5:1:11.5], 'FontSize', 8);
    set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 6.0 2.0]);
    print ('-depsc', 'results/PV_seasonal.eps');
    eps2pdf('results/PV_seasonal.eps','/usr/local/bin/gs');
end



% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
NG_price = load('traces/NG_price.csv');
%CRC = [1000/20+SR*200/4, 0.005+reshape(ones(TS,1) * NG_price',1,T)]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
CRC = [1000/20+SR*200/4, 0.005+0.10*ones(1,T)];
CRE = 443; % emmision, in gram CO2 eq

% cooling efficiency
PUE = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
GP = reshape(ones(TS,1)*load('traces/electricity_price.csv')',1,T); % grid electricity price
%GP = 0.10*ones(1,T);
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
[cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,max(plot_option + plot_demand*2),'results\tmp');
yearly_cost_grid = sum(cost_grid([2,4,5]));
total_cost_grid = yearly_cost_grid*max_payback

% net-zero, supply only
[cost_nz,capacity_nz,emission_nz] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_nz
yearly_cost_nz = sum(cost_nz([2,4,5]));
infra_cost_nz = sum(cost_nz([1,3]))*max_payback;
total_cost_nz = yearly_cost_nz*max_payback + infra_cost_nz
cost_reduction_nz = 1- total_cost_nz/total_cost_grid
[payback_nz,npv_nz] = payback_cal(infra_cost_nz,(yearly_cost_grid-yearly_cost_nz)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, integrated supply and demand opt
[cost_int,capacity_int,emission_int] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_int
yearly_cost_int = sum(cost_int([2,4,5]));
infra_cost_int = sum(cost_int([1,3]))*max_payback;
total_cost_int = yearly_cost_int*max_payback + infra_cost_int
cost_reduction_int = 1- total_cost_int/total_cost_grid
[payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% not net-zero, supply only
[cost_not,capacity_not,emission_not] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_not
yearly_cost_not = sum(cost_not([2,4,5]));
infra_cost_not = sum(cost_not([1,3]))*max_payback;
total_cost_not = yearly_cost_not*max_payback + infra_cost_not
cost_reduction_not = 1- total_cost_not/total_cost_grid
[payback_not,npv_not] = payback_cal(infra_cost_not,(yearly_cost_grid-yearly_cost_not)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
figure;
pie([cost_not(1),cost_not(3),cost_not(2),cost_not(4),cost_not(5)]*max_payback, {'PV CapEx','GE CapEx','PV OpEx','GE OpEx','Electricity'})

% not net-zero, supply + demand
[cost_notint,capacity_notint,emission_notint] = houston_grid(1,0,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
capacity_notint
yearly_cost_notint = sum(cost_notint([2,4,5]));
infra_cost_notint = sum(cost_notint([1,3]))*max_payback;
total_cost_notint = yearly_cost_notint*max_payback + infra_cost_notint
cost_reduction_notint = 1- total_cost_notint/total_cost_grid
[payback_notint,npv_notint] = payback_cal(infra_cost_notint,(yearly_cost_grid-yearly_cost_notint)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
figure;
pie([cost_notint(1),cost_notint(3),cost_notint(2),cost_notint(4),cost_notint(5)]*max_payback, {'PV CapEx','GE CapEx','PV OpEx','GE OpEx','Electricity'})

% net-zero, supply given
plot_option = 0;
payback_plot = 0;
for i = 1:1:5
    PV_capacity(i) = i*100;
    [cost_PV(i,:),capacity_PV(i,:),emission_PV(i,:)] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV_capacity(i),PV_capacity(i)],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_PV = sum(cost_PV(i,[2,4,5]));
    infra_cost_PV = sum(cost_PV(i,[1,3]))*max_payback;
    total_cost_PV(i) = yearly_cost_PV*max_payback + infra_cost_PV
    cost_reduction_PV(i) = 1- total_cost_PV(i)/total_cost_grid
    [payback_PV(i),npv_PV] = payback_cal(infra_cost_PV,(yearly_cost_grid-yearly_cost_PV)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
plot(PV_capacity,total_cost_PV);
ylabel('total expenditure ($)');
xlabel('PV capacity');
xlim([100,500]);

figure;
plot(PV_capacity,[0.1076,0.2153,0.3228,0.4116,0.4648]);
xlabel('PV capacity');
ylabel('demand powered by PV');
xlim([100,500]);

figure;
plot(PV_capacity,payback_PV);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([100,500]);

% net-zero, supply & demand
plot_option = 0;
payback_plot = 0;
for i = 1:1:5
    PV_capacity(i) = i*100;
    [cost_PVint(i,:),capacity_PVint(i,:),emission_PVint(i,:)] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV_capacity(i),PV_capacity(i)],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_PVint = sum(cost_PVint(i,[2,4,5]));
    infra_cost_PVint = sum(cost_PVint(i,[1,3]))*max_payback;
    total_cost_PVint(i) = yearly_cost_PVint*max_payback + infra_cost_PVint
    cost_reduction_PVint(i) = 1- total_cost_PVint(i)/total_cost_grid
    [payback_PVint(i),npv_PVint] = payback_cal(infra_cost_PVint,(yearly_cost_grid-yearly_cost_PVint)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
plot(PV_capacity,total_cost_PVint);
ylabel('total expenditure ($)');
xlabel('PV capacity');
xlim([100,500]);

figure;
plot(PV_capacity,[0.1754,0.3504,0.5243,0.6939,0.8285]);
xlabel('PV capacity');
ylabel('demand powered by PV');
xlim([100,500]);

figure;
plot(PV_capacity,payback_PVint);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([100,500]);
ylim([0,20]);


% net-zero, supply only, change gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:6
    gas_price(i) = i*0.02;
    [cost_gas(i,:),capacity_gas(i,:),emission_gas(i,:)] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_gas = sum(cost_gas(i,[2,4,5]));
    infra_cost_gas = sum(cost_gas(i,[1,3]))*max_payback;
    total_cost_gas(i) = yearly_cost_gas*max_payback + infra_cost_gas
    cost_reduction_gas(i) = 1- total_cost_gas(i)/total_cost_grid
    [payback_gas(i),npv_gas] = payback_cal(infra_cost_gas,(yearly_cost_grid-yearly_cost_gas)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
bar(gas_price,cost_gas*max_payback,0.2,'stacked');
ylabel('total expenditure ($)');
xlabel('gas price');
legend('PV capital','PV O&M','GE capital','GE O&M','Electricity','Location','Northwest')
xlim([0.01,0.13]);

figure;
plot(gas_price,capacity_gas(:,1),'r-',gas_price,capacity_gas(:,2),'b-',gas_price,capacity_gas(:,3),'k-');
xlabel('Gas price');
ylabel('capacity (kW)');
legend('PV','GE','Grid','Location','Northwest')
xlim([0.02,0.12]);

figure;
plot(gas_price,payback_gas);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([0.02,0.12]);
ylim([0,20]);

% net-zero, supply + demand, change gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:6
    gas_price(i) = i*0.02;
    [cost_gasint(i,:),capacity_gasint(i,:),emission_gasint(i,:)] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_gasint = sum(cost_gasint(i,[2,4,5]));
    infra_cost_gasint = sum(cost_gasint(i,[1,3]))*max_payback;
    total_cost_gasint(i) = yearly_cost_gasint*max_payback + infra_cost_gasint
    cost_reduction_gasint(i) = 1- total_cost_gasint(i)/total_cost_grid
    [payback_gasint(i),npv_gasint] = payback_cal(infra_cost_gasint,(yearly_cost_grid-yearly_cost_gasint)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
bar(gas_price,cost_gasint*max_payback,0.2,'stacked');
ylabel('total expenditure ($)');
xlabel('gas price');
legend('PV capital','PV O&M','GE capital','GE O&M','Electricity','Location','Northwest')
xlim([0.01,0.13]);

figure;
plot(gas_price,capacity_gasint(:,1),'r-',gas_price,capacity_gasint(:,2),'b-',gas_price,capacity_gasint(:,3),'k-');
xlabel('Gas price');
ylabel('capacity (kW)');
legend('PV','GE','Grid','Location','Northwest')
xlim([0.02,0.12]);

figure;
plot(gas_price,payback_gasint);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([0.02,0.12]);
ylim([0,20]);

% net-zero, supply only, gas engine 250kW, change gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:6
    gas_price(i) = i*0.02;
    [cost_gas250(i,:),capacity_gas250(i,:),emission_gas250(i,:)] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,[250,250],[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_gas250 = sum(cost_gas250(i,[2,4,5]));
    infra_cost_gas250 = sum(cost_gas250(i,[1,3]))*max_payback;
    total_cost_gas250(i) = yearly_cost_gas250*max_payback + infra_cost_gas250
    cost_reduction_gas250(i) = 1- total_cost_gas250(i)/total_cost_grid
    [payback_gas250(i),npv_gas250] = payback_cal(infra_cost_gas250,(yearly_cost_grid-yearly_cost_gas250)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
bar(gas_price,cost_gas250*max_payback,0.2,'stacked');
ylabel('total expenditure ($)');
xlabel('gas price');
legend('PV capital','PV O&M','GE capital','GE O&M','Electricity','Location','Northwest')
xlim([0.01,0.13]);

figure;
plot(gas_price,capacity_gas250(:,1),'r-',gas_price,capacity_gas250(:,2),'b-',gas_price,capacity_gas250(:,3),'k-');
xlabel('Gas price');
ylabel('capacity (kW)');
legend('PV','GE','Grid','Location','Northwest')
xlim([0.02,0.12]);

figure;
plot(gas_price,payback_gas250);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([0.02,0.12]);
ylim([0,20]);

% net-zero, supply + demand, gas engine 250kW, change gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:6
    gas_price(i) = i*0.02;
    [cost_gasint250(i,:),capacity_gasint250(i,:),emission_gasint250(i,:)] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,[250,250],[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
    yearly_cost_gasint250 = sum(cost_gasint250(i,[2,4,5]));
    infra_cost_gasint250 = sum(cost_gasint250(i,[1,3]))*max_payback;
    total_cost_gasint250(i) = yearly_cost_gasint250*max_payback + infra_cost_gasint250
    cost_reduction_gasint250(i) = 1- total_cost_gasint250(i)/total_cost_grid
    [payback_gasint250(i),npv_gasint250] = payback_cal(infra_cost_gasint250,(yearly_cost_grid-yearly_cost_gasint250)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end
figure;
bar(gas_price,cost_gasint250*max_payback,0.2,'stacked');
ylabel('total expenditure ($)');
xlabel('gas price');
legend('PV capital','PV O&M','GE capital','GE O&M','Electricity','Location','Northwest')
xlim([0.01,0.13]);

figure;
plot(gas_price,capacity_gasint250(:,1),'r-',gas_price,capacity_gasint250(:,2),'b-',gas_price,capacity_gasint250(:,3),'k-');
xlabel('Gas price');
ylabel('capacity (kW)');
legend('PV','GE','Grid','Location','Northwest')
xlim([0.02,0.12]);

figure;
plot(gas_price,payback_gasint250);
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([0.02,0.12]);
ylim([0,20]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save('results\houston_payback_seasonal.mat');
