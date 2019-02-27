close all; clear all; clc;

OFFICIAL  = 1;

addpath('functions/');

has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
% cvx_solver Gurobi;
cvx_solver SeDuMi;
cvx_precision low;
cvx_clear;

fig_path = '../figs/';
% fig_path = 'results/';
% general parameters

if ~OFFICIAL
  disp('TEST: 12 days 1 year optimization');
  N_y = 1; % Number of years
  TS = 24; % number of slots per day
  ND = 12; % number of days
else
  N_y = 7; % Number of years
  TS = 24; % number of slots per day
  ND = 12; % number of days
end
T = TS*ND; % planning period length, in hours
IP = 0.04; % idle power consumption of a server, in kW
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5*1000; % g CO2/kWh of grid power consumption
DC_power = 10000; % Data center power capacity, in kW
DC_CPU = 1000; % 1000
DC_GPU = 1000; % 1000
OP = 5/6; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
IT_CPU = DC_CPU*OP;
IT_GPU = DC_GPU*OP;
con = 0.9; % maximum utilization after consolidation
plot_PV = 0; % plot PV traces
plot_demand = 0; % plot power demand

% increasing rate
GP_increase = 1.05;%1.03; % https://www.eia.gov/forecasts/steo/report/electricity.cfm
RC_decrease = 0.85;% http://cleantechnica.com/2015/01/29/solar-costs-will-fall-40-next-2-years-heres/
CRC_increase = 1.00;% gas price:
workload_increase = 1.09; %1.09;

%%
% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.25; % batch job ratio, compared with interactive workload
IPU = 0.7/workload_increase^(N_y-1); % 0.7; % interactive workload peak utilization
au = 0.4; % average utilization of interactive workloads
IM_cpu = IT_CPU*au*IPU; % peak of interactive workload consumption, in servers
IM_gpu = IT_GPU*au*IPU; % peak of interactive workload consumption, in servers
PMR = 2.5; % peak to mean ratio, should be smaller than current PMR
%% HP Houston
% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces/SAPnew/sapTrace2.tab'};
a = zeros(IN,T);
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM_cpu(i));
end

[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties
bu = 1; % average utilization of batch job when running alone.
%% Microsoft
% interactive workloads
IN_gpu = 1; % number of interactive workloads
IL_gpu = {'traces/microsoft/unknown_cpu.csv'};

a_gpu = zeros(IN,T);
for i = 1:1:IN_gpu
    a_gpu(i,:) = interactive_process_microsoft(char(IL_gpu(i)), T, 60, 60*8, PMR, IM_gpu(i));
end
% a_gpu=(circshift(a_gpu',[12 0]))';
plot(a_gpu);

% batch job
BN_gpu= 1; % average number of batch job arrivals per timeslot
[A_gpu,BS_gpu,S_gpu,E_gpu] = batch_job_generator(T,BN_gpu*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a_gpu,2)./au)/con); % generate batch jobs based on statistical properties

%%
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
RR = [0, 1.0*DC_power]; % renewable capacity range, in kW
% RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
RC = [2150/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
% http://cleantechnica.com/2015/01/29/solar-costs-will-fall-40-next-2-years-heres/
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
CRC = [1000/20+SR*200/4, 0.005+0.06*ones(1,T)];
% CRC = [1000/20+SR*200/4, 0.005+0.07*ones(1,T)];
CRE = 443; % emmision, in gram CO2 eq

% cooling efficiency
PUE = 1.3*ones(1,T);
PUE_low = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
%GP = reshape(ones(TS,1)*load('traces\electricity_price.csv')',1,T); % grid electricity price
% GP = 0.08*ones(1,T); % grid power only as backup
GP = 0.056*ones(1,T); % grid power only as backup
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%IT
ITR = [0,DC_power*OP];
ITC = 4000/4;

GP_Array = zeros(T,N_y);

%http://www.techrepublic.com/article/state-of-the-solar-industry-10-stats-to-know/

RC_Array = zeros(1,N_y);

CRC_Array = zeros(T+1,N_y);

CRR = [0, 1400];
RR =  [0, DC_power];

% http://siliconangle.com/blog/2014/01/27/20-cloud-computing-statistics-tc0114/  
% interactive workloads & Batch job
a_Array  =  zeros(T, N_y); BS_Array =  zeros(BN*T, N_y); A_Array = repmat(A,[1 1 N_y]);

for y = 1:1:N_y    
    GP_Array(:,y) = GP*GP_increase^(y-1);
    RC_Array(1,y) = RC(1)*RC_decrease^(y-1);
    RC_Array(2,y) = RC(2);
    CRC_Array(2:T+1,y) = (CRC(2:T+1)-0.005)*CRC_increase^(y-1)+0.005;    
    CRC_Array(1,y) = CRC(1);
    
    a_Array(:, y) = a*workload_increase^(y-1);
    BS_Array(:,y)    = BS*workload_increase^(y-1);    
%     IT_C(y) = IT*workload_increase^(y-1);% Data center IT capacity, in kW
end

% https://github.com/lenhattan86/azure_trace 
% convert to GPU workload 
a_Array_gpu  =  zeros(T, N_y); BS_Array_gpu =  zeros(BN*T, N_y); A_Array_gpu = repmat(A_gpu,[1 1 N_y]);

for y = 1:1:N_y    
    a_Array_gpu(:, y) = a_gpu*workload_increase^(y-1);
    BS_Array_gpu(:,y)    = BS_gpu*workload_increase^(y-1);    
%     IT_C(y) = IT*workload_increase^(y-1);% Data center IT capacity, in kW
end

%% DR programs' parameters
% TOU: Time of Use pricing - http://cleantechnica.com/2011/12/27/time-of-day-pricing-in-texas/
%http://www.businesswire.com/news/home/20111117005294/en/TXU-Energy-Offers-Deep-Nighttime-Discounts-Electricity
night_rate = 0.05;%6.8e-2; 
peak_rate = 21.9e-2;
off_peak_rate = 0.06;%9.2e-2;
GP_TOU = off_peak_rate*ones(1,T);
for t = 1:T
    hr = mod(t,24);
    if hr <= 6 || hr >= 22 % night hours
        GP_TOU(t) = night_rate;
    elseif hr >= 13 && hr <= 18 % peak hours
        GP_TOU(t) = peak_rate;
    end
end

% Inclining Block Rates
% TODO: need to find the sources of Inclining Block Rates.
ibr_levels = [50  100]; 
ibr_rates  = [0.2 0.5];
ibr_rates = ibr_rates - mean(GP); 
ibr_len = length(ibr_levels);
for i_ibr = 2:ibr_len
    ibr_rates(i_ibr) = ibr_rates(i_ibr) - sum(ibr_rates(1:i_ibr-1));
end

% Coincident peak pricing
cpp_mth_rate = 200*mean(GP)-mean(GP);
cpp_rates = zeros(1,T);
rand_hrs = randi([11 20],24,1);
i = 1;
for t = 1:24:T
    cpp_rates(t+rand_hrs(i)-1) = cpp_mth_rate;
    i = i + 1;
end


% spinning reserve (SR) 
% Spinning Reserve From Responsive Loads http://web.ornl.gov/~webworks/cppr/y2001/rpt/116213.pdf
sr_mth_rates = [0.01 0.14 0.04 0.04 0.04  0.04  0.04 0.04 0.04 0.04 0.04 0.04];
sr_rates = zeros(1,T);
rand_hrs = randi([11 20],24,1);
i = 1;
for t = 1:24:T
    sr_rates(t+rand_hrs(i)-1) = sr_mth_rates(i);
    i = i + 1;
end

% Bilateral markets
bi_power_ratio = 1/5;
bi_rate = 0.03;

%% for prediction errors
N_samples = 50;

a_real = a; 
BS_real = BS;
R_real = R;
GP_real = GP;
CRC_real = CRC;
GP_Array_real  = GP_Array;
a_Array_real   = a_Array;
BS_Array_real  = BS_Array;
CRC_Array_real = CRC_Array;