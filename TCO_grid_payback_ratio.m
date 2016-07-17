clear all; close all;
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
init_common_settings_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
discount_factor = 0;
max_payback = 20;
payback_plot = 0;
demand_opt = 1;
scaleFactor = 1;
total_demand = (sum(sum(a))+sum(BS))/2;
a_raw = a;
BS_raw = BS;
for i = 1:1:6
    ratio(i) = (i-1)*0.2;
    BS = BS_raw/sum(BS_raw) * total_demand * ratio(i);
    a = a_raw/sum(sum(a_raw)) * total_demand * (1-ratio(i));    
    [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,0,'results\tmp');
    cost_grid = scaleFactor * cost_grid;
    capacity_grid = scaleFactor * capacity_grid;
    emission_grid = scaleFactor * emission_grid;    
    
    yearly_cost_grid = sum(cost_grid([2,4,5]));
    
    [cost(i,:),capacity(i,:),emission(i,:)] = houston_grid(demand_opt,1,demand_opt,demand_opt,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
    cost(i,:) = scaleFactor * cost(i,:);
    capacity(i,:) = scaleFactor * capacity(i,:);
    emission(i,:) = scaleFactor * emission(i,:); 
    
    yearly_cost = sum(cost(i,[2,4,5]));
    infra_cost = sum(cost(i,[1,3]))*max_payback;    
    [payback(i),npv_nz(i)] = payback_cal(infra_cost,(yearly_cost_grid-yearly_cost)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
bar(ratio,cost,0.2,'stacked');
ylabel('annual expenditure ($)');
xlabel('flexible workload ratio');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill','Location','Northeast');
ylim([0,max(sum(cost,2))*1.3]);
xlim([-0.1,1.1]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', [fig_path 'cost_batch.eps']);
% saveas(gcf,'results\expenditure_comp_pb_gp.fig');

figure;
plot(ratio,capacity(:,1),'r-',ratio,capacity(:,2),'b-',ratio,capacity(:,3),'k-');
xlabel('flexible workload ratio');
ylabel('capacity (kW)');
legend('PV','GE','Grid');
ylim([0,max(max(capacity))*1.3]);
xlim([0,1.0]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', [fig_path 'capacity_batch.eps']);
% saveas(gcf,'results\capacity_comp_pb_gp.fig');

% figure;
% plot(ratio,payback)
% xlabel('flexible workload ratio');
% ylabel('pay back period (year)');
% ylim([0,max_payback]);
% xlim([0,1.0]);

save('results\houston_payback_ratio.mat');

%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}