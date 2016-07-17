close all; clear all;

addpath('functions/');

has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_solver Gurobi;
cvx_precision low;
cvx_clear;

% fig_path = '../figs/';
fig_path = 'results/';
% general parameters
TS = 24; % number of slots per day
ND = 12; % number of days
T = TS*ND; % planning period length, in hours
IP = 0.04; % idle power consumption of a server, in kW
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5*1000; % g CO2/kWh of grid power consumption
DC_power = 1000; % Data center power capacity, in kW
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
RR = [0, 1.0*DC_power]; % renewable capacity range, in kW
% RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
RC = [2150/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
% http://cleantechnica.com/2015/01/29/solar-costs-will-fall-40-next-2-years-heres/
RE = 32; % emmision, in gram CO2 eq for renewable

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
end

% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
NG_price = load('traces/NG_price.csv');
%CRC = [1000/20+SR*200/4, 0.005+reshape(ones(TS,1) * NG_price',1,T)]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
CRC = [1000/20+SR*200/4, 0.005+0.06*ones(1,T)];
% CRC = [1000/20+SR*200/4, 0.005+0.07*ones(1,T)];
CRE = 443; % emmision, in gram CO2 eq % gas emission.

% cooling efficiency
PUE = 1.3*ones(1,T);
PUE_low = 1.2*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
%GP = reshape(ones(TS,1)*load('traces\electricity_price.csv')',1,T); % grid electricity price
% GP = 0.08*ones(1,T); % grid power only as backup
GP = 0.056*ones(1,T); % grid power only as backup
BP = 0.00*ones(1,T); % sell back price
RP = 0.00*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%IT
ITR = [0,DC_power*OP];
ITC = 4000/4;

%% DR programs' parameters
% TOU: Time of Use pricing - http://cleantechnica.com/2011/12/27/time-of-day-pricing-in-texas/
%http://www.businesswire.com/news/home/20111117005294/en/TXU-Energy-Offers-Deep-Nighttime-Discounts-Electricity
% night_rate = 0; 
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
% cpp_mth_rate = 200*mean(GP);
% cpp_mth_rate = cpp_mth_rate/30;
cpp_mth_rate = 200*mean(GP);
cpp_mth_rate = cpp_mth_rate/4;
cpp_rates = zeros(1,T);
rand_hrs = randi([11 20],24,1);
% rand_hrs = randi([16 17],T/TS,1);
rand_hrs(1) = 16;
i = 1;
for t = 1:24:24
    cpp_rates(t+rand_hrs(i)-1) = cpp_mth_rate;
    i = i + 1;
end


% spinning reserve (SR) 
% Spinning Reserve From Responsive Loads http://web.ornl.gov/~webworks/cppr/y2001/rpt/116213.pdf
sr_mth_rates = [0.01 0.14 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04];
sr_rates = zeros(1,T);
rand_hrs = randi([12 18],T/TS,1);
i = 1;
rand_hrs(1) = 17;
for t = 1:24:T
    sr_rates(t+rand_hrs(i)-1) = sr_mth_rates(i);
    i = i + 1;
end

% Bilateral markets
bi_power_ratio = 1/5;
bi_rate = 0.05;