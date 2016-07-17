% based on TCO_grid_electricity_0115.m
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
load_exist = 0;
mat_location = 'results/TCO_grid_electricity_0115.mat';

    
%% Load the common settings.
init_common_settings;

%% additional parameters by Tan
NY = 7; % number of Years.
GP_increase = 1.05; % https://www.eia.gov/forecasts/steo/report/electricity.cfm
GP_Array = zeros(NY,1);

%http://www.techrepublic.com/article/state-of-the-solar-industry-10-stats-to-know/
RC_decrease = 1;% http://cleantechnica.com/2015/01/29/solar-costs-will-fall-40-next-2-years-heres/
RC_Array = zeros(NY,1);

CRC_increase = 1.00;% http://blog.rmi.org/blog_2014_08_28_helping_battery_cost_declines_keep_going_and_going
CRC_Array = zeros(NY,T);

CRR = [0,1400];
RR = [0, 2*DC_power];

% http://siliconangle.com/blog/2014/01/27/20-cloud-computing-statistics-tc0114/  
% interactive workloads & Batch job
workload_increase = 1.09;
a_Array  =  zeros(NY, IN,T);
BS_Array =  zeros(NY, BN*T);

for y = 1:1:NY    
    GP_Array(y) = mean(GP)*GP_increase^(y-1);
    RC_Array(y) = RC(1)*RC_decrease^(y-1);
    CRC_Array(y,:) = (CRC(2:T+1)-0.005)*CRC_increase^(y-1)+0.005;    

%     aTemp = interactive_process(char(IL(i))^(y-1), T, 12, 4, PMR, IM(i));
    a_Array(y, :, :) = a*workload_increase^(y-1);
    BS_Array(y,:)    = BS*workload_increase^(y-1);
    IT_C(y) = IT*workload_increase^(y-1);% Data center IT capacity, in kW
%     IT_C(y) = 999999;% Data center IT capacity, in kW
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
discount_factor = 0;
demand_opt = 1;
nz = 0;

capacity_prev = [0 0 0];
for y = 1:1:NY
    GP = GP_Array(y)*ones(1,T);
    RC(1) = RC_Array(y);
    CRC(2:T+1) = CRC_Array(y,:);
    a = squeeze(a_Array(y,:,:))';
    BS = squeeze(BS_Array(y,:))';   
    [cost(y,:),capacity(y,:),emission(y,:)] = houston_grid_annual(demand_opt,nz,demand_opt, ...
        demand_opt,T,IP,PP,OP,IT_C(y),CO2_grid,a,au,con,BS,A,S,E,bu, ...
        PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results/temp', capacity_prev);  
    installed_capacity(y,:) = capacity(y,:) - capacity_prev;
%     installed_capacity(y,3) = capacity(y,3);
    capacity_prev = capacity(y,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all;
figure;
bar(cost,0.2,'stacked');
ylabel('annual expenditure ($)');
xlabel('year');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill');
ylim([0,sum(cost(1,:))*1.5]);
xlim([0,NY+1]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 3.2]);
print ('-depsc', [fig_path 'annual_change_cost.eps']);

figure;
bar(installed_capacity);
ylabel('capacity (kW)');
xlabel('year');
legend('PV','GE','Grid','Location','northwest');
ylim([min(min(installed_capacity))*1.05,max(max(installed_capacity))*1.3]);
xlim([0,NY+1]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 3.2]);
print ('-depsc', [fig_path 'annual_change_capacity.eps']);

save('results/TCO_annual_change.mat');