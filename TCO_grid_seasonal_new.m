close all;
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
DC_power = 300; % Data center power capacity, in kW
OP = 5/6; % over-provisioning ratio of IT and power
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
R_raw = zeros(12,24*31);
PVN = zeros(12,1);
CF = zeros(12,1);
RCF = load('traces/PV_CF.csv');
for i = 1:1:12
    tmp = load(strcat('traces/Houston/',int2str(i),'.csv'));
    PVN(i) = size(tmp,1)/12;
    for t = 1:1:PVN(i)
        R_raw(i,t) = max(0,mean(tmp((t-1)*12+1:(t-1)*12+12,6)));
    end
end
for i = 1:1:12
    CF(i) = mean(R_raw(i,:))/max(max(R_raw));
end
for i = 1:1:12
    R((i-1)*24+1:(i-1)*24+24) = sum(reshape(R_raw(i,1:PVN(i)),24,PVN(i)/24),2);
end
R = R/max(R);
RR = [0, 500]; % renewable capacity range, in kW
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
CRR = [250,250]; % controled renewable capacity range, in kW
NG_price = load('traces/NG_price.csv');
%CRC = [1000/20+SR*200/4, 0.005+reshape(ones(TS,1) * NG_price',1,T)]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
CRC = [1000/20+SR*200/4, 0.005+0.10*ones(1,T)];
CRE = 443; % emmision, in gram CO2 eq

% cooling efficiency
PUE = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
%GP = reshape(ones(TS,1)*load('traces\electricity_price.csv')',1,T); % grid electricity price
GP = 10*ones(1,T); % grid power only as backup
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

% net-zero, supply only, change pv capacity and gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:3
    gas_price(i) = i*0.02+0.03;
    for j = 1:1:10
        PV(j) = j*50;
        [cost_gasPV(i,j,:),capacity_gasPV(i,j,:),emission_gasPV(i,j,:)] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV(j),PV(j)],RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
        yearly_cost_gasPV = sum(cost_gasPV(i,j,[2,4,5]));
        infra_cost_gasPV = sum(cost_gasPV(i,j,[1,3]))*max_payback;
        total_cost_gasPV(i,j) = yearly_cost_gasPV*max_payback + infra_cost_gasPV
    end
end

figure;
for i = 1:1:3
    hold all
    plot(PV,total_cost_gasPV(i,:));
end
box;
ylabel('total expenditure ($)');
xlabel('PV capacity (kW)');
legend('gas price = $0.05/kWh','gas price = $0.07/kWh','gas price = $0.09/kWh','Location','Northwest')
xlim([50,500]);

% net-zero, supply & demand, change pv capacity and gas price
plot_option = 0;
payback_plot = 0;
for i = 1:1:3
    gas_price(i) = i*0.02+0.03;
    for j = 1:1:10
        PV(j) = j*50;
        [cost_gasPVint(i,j,:),capacity_gasPVint(i,j,:),emission_gasPVint(i,j,:)] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV(j),PV(j)],RC,RE,SR,CRR,[1000/20+SR*200/4, 0.005+gas_price(i)*ones(1,T)],CRE,P,GP,RP,BP,plot_option,'results\tmp');
        yearly_cost_gasPVint = sum(cost_gasPVint(i,j,[2,4,5]));
        infra_cost_gasPVint = sum(cost_gasPVint(i,j,[1,3]))*max_payback;
        total_cost_gasPVint(i,j) = yearly_cost_gasPVint*max_payback + infra_cost_gasPVint
    end
end

figure;
for i = 1:1:3
    hold all
    plot(PV,total_cost_gasPVint(i,:));
end
box;
ylabel('total expenditure ($)');
xlabel('PV capacity (kW)');
legend('gas price = $0.05/kWh','gas price = $0.07/kWh','gas price = $0.09/kWh','Location','Northwest')
xlim([0.01,0.13]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save('results\houston_payback_seasonal.mat');
