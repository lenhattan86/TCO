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

gs_path = '/usr/local/bin/gs'; % path of pdf files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_exist = 0;
mat_location = 'results/houston_payback_seasonal_0213.mat';

if load_exist == 1
    load('results/houston_payback_seasonal_0213.mat');
else
    %% Load the common settings.
    init_common_settings;
    GP = 0.07*ones(1,T); %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          SOLVER                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_option = 1;
    discount_factor = 0;
    max_payback = 20;
    payback_plot = 0;
    nz = 0;
    
    % grid only
    [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,[fig_path 'nz-grid']);
    yearly_cost_grid = sum(cost_grid([2,4,5]));
    total_cost_grid = yearly_cost_grid*max_payback
    
    % net-zero, supply only
    [cost_supply,capacity_supply,emission_supply] = houston_grid(0,nz,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,[fig_path 'nz-supply']);
    yearly_cost_supply = sum(cost_supply([2,4,5]));
    infra_cost_supply = sum(cost_supply([1,3]))*max_payback;
    total_cost_supply = yearly_cost_supply*max_payback + infra_cost_supply
    cost_reduction_supply = 1- total_cost_supply/total_cost_supply
    [payback_supply,npv_supply] = payback_cal(infra_cost_supply,(yearly_cost_grid-yearly_cost_supply)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
    
    % net-zero, demand only
    [cost_demand,capacity_demand,emission_demand] = houston_grid(1,nz,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE_low,R(1,:),0.7*[DC_power,DC_power],RC,RE,SR,0.5*[DC_power,DC_power],CRC,CRE,P,GP,RP,BP,plot_option,[fig_path 'nz-demand']);
    yearly_cost_demand = sum(cost_demand([2,4,5]));
    infra_cost_demand = sum(cost_demand([1,3]))*max_payback;
    total_cost_demand = yearly_cost_demand*max_payback + infra_cost_demand
    cost_reduction_demand = 1- total_cost_demand/total_cost_demand
    [payback_demand,npv_demand] = payback_cal(infra_cost_demand,(yearly_cost_grid-yearly_cost_demand)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
    
    % net-zero, integrated
    [cost_int,capacity_int,emission_int] = houston_grid(1,nz,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE_low,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,[fig_path 'nz-int']);
    yearly_cost_int = sum(cost_int([2,4,5]));
    infra_cost_int = sum(cost_int([1,3]))*max_payback;
    total_cost_int = yearly_cost_int*max_payback + infra_cost_int
    cost_reduction_int = 1- total_cost_int/total_cost_int
    [payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
    
end

save('results/houston_payback_seasonal_0213.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
return; 

figure;
bar([cost_grid;cost_supply;cost_demand;cost_int],0.2,'stacked');
ylabel('annual expenditure ($)');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill');
ylim([0,sum(cost_grid)*1.3]);
set(gca,'xticklabel',{'Grid','Supply','Demand','Integrated'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
print ('-depsc', [fig_path 'cost_comp.eps']);
% eps2pdf('results/cost_comp.eps',gs_path);

figure;
bar([capacity_grid;capacity_supply;capacity_demand;capacity_int]);
ylabel('capacity (kW)');
legend('PV','GE');
ylim([0,max(capacity_demand)*1.3]);
set(gca,'xticklabel',{'Grid only','Supply opt','Supply given','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
print ('-depsc', [fig_path 'capacity_comp.eps']);
% eps2pdf('results/capacity_comp.eps',gs_path);

figure;
bar([emission_grid;emission_supply;emission_demand;emission_int]/1000,0.2,'stacked');
ylabel('annual emission (kg)');
% ylim([0,max(emission_grid)*1.3]);
set(gca,'xticklabel',{'Grid only','Supply opt','Supply given','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
print ('-depsc', [fig_path 'emission_comp.eps']);
% eps2pdf('results/emission_comp.eps',gs_path);

