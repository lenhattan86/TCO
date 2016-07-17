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

load_exist = 0;
mat_location = 'results/TCO_grid_IT_0115.mat';

if load_exist == 1
    load('results/TCO_grid_IT_0115.mat');
else

init_common_settings;

%IT
ITR = [0,DC_power*OP];
ITC = 4000/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 1;
discount_factor = 0;
max_payback = 20;
payback_plot = 0;
nz = 0;
% grid only
RC = [2600/20, 0.005+0];RR=[0,500];CRC = [1000/20+SR*200/4, 0.005+0.06];GP = 0.08*ones(1,T);CRR=[0,1400];
[cost_grid,capacity_grid,emission_grid] = houston_grid_IT(0,0,0,0,T,IP,PP,DC_power,[DC_power*OP,DC_power*OP],ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,plot_option,'results/tmp');
yearly_cost_grid = sum(cost_grid([2,4,5]));

% net-zero, supply only
[cost_nz,capacity_nz,emission_nz] = houston_grid_IT(0,nz,0,0,T,IP,PP,DC_power,ITR,ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results/tmp');
yearly_cost_nz = sum(cost_nz([2,4,5]));
infra_cost_nz = sum(cost_nz([1,3]))*max_payback;
[payback_nz,npv_nz] = payback_cal(infra_cost_nz,(yearly_cost_grid-yearly_cost_nz)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, demand only
[cost_ds,capacity_ds,emission_ds] = houston_grid_IT(1,nz,1,1,T,IP,PP,DC_power,[DC_power*OP,DC_power*OP],ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[150,150],RC,RE,SR,[633,633],CRC,CRE,P,GP,RP,BP,plot_option,'results/tmp');
yearly_cost_ds = sum(cost_ds([2,4,5]));
infra_cost_ds = sum(cost_ds([1,3]))*max_payback;
[payback_ds,npv_ds] = payback_cal(infra_cost_ds,(yearly_cost_grid-yearly_cost_ds)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

% net-zero, integrated supply and demand opt
[cost_int,capacity_int,emission_int] = houston_grid_IT(1,nz,1,1,T,IP,PP,DC_power,ITR,ITC,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),RR,RC,RE,SR,CRR,CRC,CRE,P,GP,RP,BP,plot_option,'results/tmp');
yearly_cost_int = sum(cost_int([2,4,5]));
infra_cost_int = sum(cost_int([1,3]))*max_payback;
[payback_int,npv_int] = payback_cal(infra_cost_int,(yearly_cost_grid-yearly_cost_int)*ones(1,max_payback),discount_factor,max_payback,payback_plot)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;
bar([cost_grid;cost_nz;cost_ds;cost_int],0.2,'stacked')
ylim([0,sum(cost_grid)*2]);
ylabel('annual expenditure ($)');
legend('PV infrastructure cost','PV O&M cost','GE infrastructure cost','GE O&M cost','Electricity bill','IT capital cost','Location','Northeast');
set(gca,'xticklabel',{'Grid','Supply','Demand','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cost_IT.eps');
eps2pdf('results/cost_IT.eps','/usr/local/bin/gs');

figure;
bar([capacity_grid;capacity_nz;capacity_ds;capacity_int])
legend('PV','GE','Grid','IT');
ylabel('capacity (kW)');
set(gca,'xticklabel',{'Grid','Supply','Demand','Our solution'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/capacity_IT.eps');
eps2pdf('results/capacity_IT.eps','/usr/local/bin/gs');
%{
figure;
bar([cost_nz(1),cost_nz(3),0,cost_nz(6);cost_nz(2),cost_nz(4),cost_nz(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_ds(1),cost_ds(3),0,cost_ds(6);cost_ds(2),cost_ds(4),cost_ds(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');

figure;
bar([cost_int(1),cost_int(3),0,cost_int(6);cost_int(2),cost_int(4),cost_int(5),0]*max_payback)
legend('PV','GE','Grid','IT');
%ylim([0,sum(cost_nz)*max_payback]);
set(gca,'xticklabel',{'CapEx','OpEx'});
ylabel('total cost ($)');
%}
save('results/TCO_grid_IT_0115.mat');