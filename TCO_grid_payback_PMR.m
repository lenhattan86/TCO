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
init_common_settings_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
discount_factor = 0;
max_payback = 20;
payback_plot = 0;
demand_opt = 0;
for i = 1:1:6
    PMR(i) = 1+i*0.3; % peak to mean ratio, should be smaller than current PMR
    for j = 1:1:IN
        a(j,:) = interactive_process(char(IL(j)), T, 12, 4, PMR(i), IM(j)*PMR(i)/3);
    end
    workload(i) = sum(sum(a));
%     CRC = [1000/20+SR*200/4, 0.005+0.06];GP = 0.08*ones(1,T);RC = [2600/20, 0.005+0];    
%     RR=[0,1500];
%     CRR=[0,1400];
    [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid, ... 
        a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,0,'results\tmp');
    yearly_cost_grid = sum(cost_grid([2,4,5]));
%     [cost(i,:),capacity(i,:),emission(i,:)] = houston_grid(demand_opt,1,demand_opt,demand_opt,...
%         T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,...
%         BP,plot_option,'results\tmp');
    [cost(i,:),capacity(i,:),emission(i,:)] = houston_grid(demand_opt,1,demand_opt,demand_opt, ... 
        T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,...
        BP,plot_option,'results\tmp');
    
    yearly_cost = sum(cost(i,[2,4,5]));
    infra_cost = sum(cost(i,[1,3]))*max_payback;
    [payback(i),npv_nz(i)] = payback_cal(infra_cost,(yearly_cost_grid-yearly_cost)*ones(1,max_payback),discount_factor,max_payback,payback_plot)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all;
figure;
bar(PMR,cost,0.2,'stacked');
ylabel('annual expenditure ($)');
xlabel('Peak-to-Mean Ratio');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill','Location','Northwest');
ylim([0,max(sum(cost,2))*1.3]);
% xlim([0.1,0.32]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', [fig_path 'cost_PMR.eps']);
% saveas(gcf,'results\expenditure_comp_pb_gp.fig');

figure;
plot(PMR,capacity(:,1),'r-',PMR,capacity(:,2),'b-',PMR,capacity(:,3),'k-');
xlabel('Peak-to-Mean Ratio');
ylabel('capacity (kW)');
legend('PV','GE','Grid');
ylim([0,max(max(capacity))*1.3]);
% xlim([0.12,0.3]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', [fig_path 'capacity_PMR.eps']);
% saveas(gcf,'results\capacity_comp_pb_gp.fig');

% figure;
% plot(0.12:0.02:0.3,payback)
% xlabel('Peak-to-Mean Ratio');
% ylabel('pay back period (year)');
% ylim([0,max_payback]);
% xlim([0.12,0.3]);

save('results\houston_payback_PMR.mat');

%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}