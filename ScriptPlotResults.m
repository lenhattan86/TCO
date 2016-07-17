

%% Plot results

%% emissions per country vs. data centers
figure_settings;
yArray = [80 142 146 178];
bar(yArray, 0.2);
ylabel('mega tons CO_2');
%xlabel('countries','FontSize',fontAxis);
xLabels = {'data centers', 'Argentina', 'Netherlands', 'Malaysia'};
%legend(legendStr,'Location','southeast','FontSize',fontLegend);


%% Prediction errors TCO_comparison_err_x.xx
figure_settings;
fileList = [ 'results/TCO_comparison_err_0.mat   ' ...
            ;'results/TCO_comparison_err_0.05.mat' ...
            ;'results/TCO_comparison_err_0.1.mat ' ...
            ;'results/TCO_comparison_err_0.15.mat' ...
            ;'results/TCO_comparison_err_0.2.mat ' ...
            ;'results/TCO_comparison_err_0.25.mat' ...
            ;'results/TCO_comparison_err_0.3.mat ' ...
            ;'results/TCO_comparison_err_0.35.mat' ...
            ;'results/TCO_comparison_err_0.4.mat ' ...
            ;'results/TCO_comparison_err_0.45.mat' ...
            ;'results/TCO_comparison_err_0.5.mat '...
            ;'results/TCO_comparison_err_0.55.mat'...
            ;'results/TCO_comparison_err_0.6.mat '...
            ;'results/TCO_comparison_err_0.65.mat'...
            ;'results/TCO_comparison_err_0.7.mat '...
            ];
        
fileLen = length(fileList(:,1));
iFile = 0;
for iCount = 1:fileLen
    if exist(fileList(iCount,:),'file')
        iFile = iFile + 1;
        load(fileList(iCount,:));
%         predErrStd
        normalizedErrors(iFile) = predErrStd; 
%         capacity_supply_mean
        cost_sum_supply(iFile) = sum(sum(cost_supply_mean));
        cost_sum_supply_actual(iFile) = sum(sum(cost_supply_mean_actual));
        emission_sum_supply_actual(iFile) = sum(sum(emission_supply_mean_actual));
        
        cost_sum_int(iFile) = sum(sum(cost_int_mean));
        cost_sum_int_actual(iFile) = sum(sum(cost_int_mean_actual));
        emission_sum_int_actual(iFile) = sum(sum(emission_int_mean_actual));
%         capacity_int_mean
    end
end
   
%figure;

%ylabel('capacity (kW)','FontSize',fontAxis);

% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_sum_int)*1.1]);
% % xlim([0.5 4.5]);
% set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [fig_path 'capacity_err_comparison.eps']);

load('results/TCO_comparison.mat');

cost_sum_grid   = sum(sum(cost_grid'));
cost_sum_demand = sum(sum(cost_demand'));

emission_sum_grid = sum(sum(emission_grid'));
emission_sum_demand = sum(sum(emission_demand'));

% figure;
% bar([capacity_supply_mean;capacity_sum_supply;capacity_sum_demand;capacity_sum_int],1);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_sum_int)*1.1]);
% xlim([0.5 4.5]);
% set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [fig_path 'capacity_err_comparison.eps']);

figure;
plot(normalizedErrors, cost_sum_grid*ones(iFile,1), gridLine, 'linewidth', lineWidth);
hold on;
plot(normalizedErrors, cost_sum_supply_actual, supLine, 'linewidth', lineWidth);
hold on;
plot(normalizedErrors, cost_sum_demand*ones(iFile,1), demLine, 'linewidth', lineWidth);
hold on;
plot(normalizedErrors, cost_sum_int_actual, propLine, 'linewidth', lineWidth);
ylabel('expenditure ($)');
xlabel('normalized errors','FontSize',fontAxis);
ylim([0,max(cost_sum_grid)*1.1]);
xlim([0,max(normalizedErrors)]);
% legendStr = {gridStr,demandStr};
legendStr = {gridStr,supplyStr,demandStr,integratedStr};

% ylim([0,max(cost_sum_demand)*1.1]);
% legendStr = {demandStr,integratedStr};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_err_comparison.eps']);

figure;
plot(normalizedErrors, emission_sum_grid*ones(iFile,1)/1000,gridLine,  'linewidth', lineWidth);
hold on;
plot(normalizedErrors, emission_sum_supply_actual/1000, supLine, 'linewidth', lineWidth);
hold on;
plot(normalizedErrors, emission_sum_demand*ones(iFile,1)/1000, demLine, 'linewidth', lineWidth);
hold on;
plot(normalizedErrors, emission_sum_int_actual/1000, propLine, 'linewidth', lineWidth);
ylabel('emissions (kg)');
xlabel('normalized errors','FontSize',fontAxis);
ylim([0,max(emission_sum_grid/1000)*1.1]);
xlim([0,max(normalizedErrors)]);
% xlim([0.5 4.5]);
% legendStr = {gridStr,demandStr};
legendStr = {gridStr,supplyStr,demandStr,integratedStr};
% legendStr = {demandStr,integratedStr};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'emission_err_comparison.eps']);

PROP_costSaving = cost_sum_int_actual(length(cost_sum_int_actual))/cost_sum_grid;
PROP_costSaving

PROP_emissionReduction = emission_sum_int_actual(length(emission_sum_int_actual))/emission_sum_grid;
PROP_emissionReduction


%% long-term TCO_comparison
figure_settings;
load('results/TCO_comparison.mat');
capacity_sum_grid   = max(capacity_grid');
capacity_sum_supply = max(capacity_supply');
capacity_sum_demand = max(capacity_demand');
capacity_sum_int    = max(capacity_int');

cost_sum_grid   = sum(cost_grid');
cost_sum_supply = sum(cost_supply');
cost_sum_demand = sum(cost_demand');
cost_sum_int    = sum(cost_int');

emission_sum_grid = sum(emission_grid');
emission_sum_supply = sum(emission_supply');
emission_sum_demand = sum(emission_demand');
emission_sum_int = sum(emission_int');
   
figure;
bar([capacity_sum_grid;capacity_sum_supply;capacity_sum_demand;capacity_sum_int],1);
ylabel('capacity (kW)','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(capacity_sum_int)*1.1]);
xlim([0.5 4.5]);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_comparison.eps']);

figure;
dumpCosts = zeros(1,5);
h  = bar([cost_sum_grid;cost_sum_supply;cost_sum_demand;cost_sum_int],barWidth*2/3,'stacked');
%h  = bar([dumpCosts;dumpCosts;dumpCosts;dumpCosts],barWidth*2/3,'stacked');
ylabel('expenditure ($)');
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
xlim([0.5 5]);
ylim([0,sum(cost_sum_grid)*1.1]);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
if IS_CHANGE_COLOR
    for hi = 1:numel(h)
        set(h(hi), 'FaceColor', colorset(hi,:))
    end
end
print ('-depsc', [FIG_PATH 'cost_comparison.eps']);

dumpEmissions = zeros(1,5);
figure;
bar([emission_sum_grid;emission_sum_supply;emission_sum_demand;emission_sum_int]/1000*0,barWidth*2/3,'stacked');
%bar([dumpEmissions;dumpEmissions;dumpEmissions;dumpEmissions]/1000,barWidth*2/3,'stacked');
ylabel('emissions (kg)');
xlim([0.5 4.5]);
ylim([0,sum(emission_sum_grid/1000)*1.1]);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'emission_comparison.eps']);

%% TCO_long_term
figure_settings;
load('results/TCO_long_term.mat');
   
figure;
bar(capacity',1);
ylabel('capacity (kW)');
xlabel('year','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([min(min(capacity))*1.05,max(max(capacity))*1.05]);
xlim([0.5,N_y+0.5]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'long_term_capacity.eps']);

cost = cost';
figure;
bar(cost,barWidth,'stacked');
ylabel('annual expenditure ($)');
xlabel('year','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
ylim([0,sum(cost(N_y,:))*1.05]);
xlim([0.5,N_y+0.5]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'long_term_cost.eps']);

emissions = emission'/1000;
figure;
bar(emissions,barWidth,'stacked');
ylabel('annual emissions (kg)');
xlabel('year','FontSize',fontAxis);
xlim([0.5,N_y+0.5]);
ylim([0,sum(emissions(1,:))*1.05]);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [fig_path 'long_term_co2.eps']);

norm_emission = sum(emissions,2)./total_demand';
figure;
bar(norm_emission,barWidth);
ylabel('normalized emissions (g/kWh)');
xlabel('year','FontSize',fontAxis);
xlim([0.5,N_y+0.5]);
ylim([0,max(norm_emission)*1.05]);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'long_term_co2_norm.eps']);

%% Annual change
figure_settings;
load('results/TCO_annual_change.mat');
   
figure;
bar(installed_capacity,1);
ylabel('capacity (kW)');
xlabel('year','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','north','FontSize',fontLegend);
ylim([min(min(installed_capacity))*1.05,max(max(installed_capacity))*1.05]);
xlim([0.5,NY+0.5]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'annual_change_capacity.eps']);

figure;
bar(cost,barWidth,'stacked');
ylabel('annual expenditure ($)');
xlabel('year','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
ylim([0,sum(cost(NY,:))*1.05]);
xlim([0.5,NY+0.5]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'annual_change_cost.eps']);

figure;
bar(emission/1000,barWidth,'stacked');
ylabel('annual emissions (kg)');
xlabel('year','FontSize',fontAxis);
xlim([0.5,NY+0.5]);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,max(emission_grid)*1.3]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'annual_change_co2.eps']);

%%
figure_settings;
load('results/houston_payback_seasonal_0213.mat');

figure;
bar([capacity_grid;capacity_supply;capacity_demand;capacity_int],1);
ylabel('capacity (kW)','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(capacity_int)*1.1]);
xlim([0.5 4.5]);
set(gca,'xticklabel',{gridOnlyStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [fig_path 'capacity_comp.eps']);
% eps2pdf('results/capacity_comp.eps',gs_path);

figure;
bar([cost_grid;cost_supply;cost_demand;cost_int],barWidth*2/3,'stacked');
ylabel('annual expenditure ($)');
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
xlim([0.5 4.5]);
ylim([0,sum(cost_grid)*1.05]);
set(gca,'xticklabel',{gridOnlyStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [fig_path 'cost_comp.eps']);

figure;
bar([emission_grid;emission_supply;emission_demand;emission_int]/1000,barWidth*2/3,'stacked');
ylabel('annual emissions (kg)');
xlim([0.5 4.5]);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,max(emission_grid)*1.3]);
set(gca,'xticklabel',{gridOnlyStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [fig_path 'emission_comp.eps']);
% eps2pdf('results/emission_comp.eps',gs_path);


%% IT capacity
figure_settings;
load('results/TCO_grid_IT_0115.mat');

figure;
bar([capacity_grid;capacity_nz;capacity_ds;capacity_int],1)
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr, ITCapacityStr};
legend(legendStr,'Location','north','FontSize',fontLegend);
ylabel('capacity (kW)');
ylim([0 max(capacity_grid)*1.4]);
xlim([0.5 4.5]);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_IT.eps']);

figure;
bar([cost_grid;cost_nz;cost_ds;cost_int],barWidth*2/3,'stacked');
ylim([0,sum(cost_grid)*1.05]);
xlim([0.5 4.5]);
ylabel('annual expenditure ($)','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','west','FontSize',fontLegend);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_IT.eps']);

%% Demand Response
figure_settings;
load('results/TCO_DR_evaluation.mat');

%Without DR

% figure;
% bar([capacity_supply;capacity_demand;capacity_int],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_NoDR.eps']);
% 
% figure;
% bar([cost_supply;cost_demand;cost_int],barWidth,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply)*2]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_NoDR.eps']);
% 
% 
% %Time of Use
% figure;
% bar([cost_supply_tou;cost_demand_tou;cost_int_tou],barWidth/2,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply_tou)*2]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_tou.eps']);
% 
% figure;
% bar([capacity_supply_tou;capacity_demand_tou;capacity_int_tou],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand_tou)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_tou.eps']);
% 
% %CPP
% figure;
% bar([cost_supply_cpp;cost_demand_cpp;cost_int_cpp],barWidth/2,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply_cpp)*2]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_cpp.eps']);
% 
% figure;
% bar([capacity_supply_cpp;capacity_demand_cpp;capacity_int_cpp],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand_cpp)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_cpp.eps']);
% 
% % SR
% figure;
% bar([cost_supply_sr;cost_demand_sr;cost_int_sr],barWidth/2,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply_sr)*2]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_sr.eps']);
% 
% figure;
% bar([capacity_supply_sr;capacity_demand_sr;capacity_int_sr],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand_sr)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_sr.eps']);
% 
% %ibr
% figure;
% bar([cost_supply_ibr;cost_demand_ibr;cost_int_ibr],barWidth/2,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply_ibr)*2]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_ibr.eps']);
% 
% figure;
% bar([capacity_supply_ibr;capacity_demand_ibr;capacity_int_ibr],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand_ibr)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_ibr.eps']);
% 
% %Wholesale
% figure;
% bar([cost_supply_BiMarkets;cost_demand_BiMarkets;cost_int_BiMarkets],barWidth/2,'stacked');
% ylabel('annual expenditure ($)','FontSize',fontAxis);
% legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
% legend(legendStr,'Location','northeast','FontSize',fontLegend);
% ylim([0,sum(cost_supply_BiMarkets)*2.5]);
% xlim([0.5,3.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'cost_comp_BiMarkets.eps']);
% 
% figure;
% bar([capacity_supply_BiMarkets;capacity_demand_BiMarkets;capacity_int_BiMarkets],groupBarWidth);
% ylabel('capacity (kW)','FontSize',fontAxis);
% legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
% legend(legendStr,'Location','northwest','FontSize',fontLegend);
% ylim([0,max(capacity_demand)*1.5]);
% set(gca,'xticklabel',{supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
% set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
% print ('-depsc', [FIG_PATH 'capacity_comp_BiMarkets.eps']);

%% Plot the comparisons of Capacity, Cost, Emission in DR programs
figure_settings;
load('results/TCO_DR_evaluation.mat');

% % 1. Supply only
% % 1.a. Capacity
%     figure;
%     bar([capacity_supply; capacity_supply_tou; capacity_supply_ibr; ...
%         capacity_supply_cpp; capacity_supply_sr; capacity_supply_BiMarkets],groupBarWidth);
%     ylabel('capacity (kW)','FontSize',fontAxis);    
%     legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
%     legend(legendStr,'Location','northwest','FontSize',fontLegend);
%     ylim([0,max(capacity_demand)*1.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'capacity_DR_sup_only.eps']);    
%     
% % 1.b. Cost
%     figure;
%     bar([cost_supply;cost_supply_tou;cost_supply_ibr;cost_supply_cpp; ... 
%         cost_supply_sr; cost_supply_BiMarkets],barWidth,'stacked');
%     ylabel('annual expenditure ($)','FontSize',fontAxis);
%     legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
%     legend(legendStr,'Location','northeast','FontSize',fontLegend);
%     ylim([0,sum(cost_supply)*2.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'cost_DR_sup_only.eps']);
% % 1.c. Emission
%     figure;
%     bar([emission_supply;emission_supply_tou;emission_supply_ibr; ...
%         emission_supply_cpp;emission_supply_sr; emission_supply_BiMarkets ],barWidth,'stacked');
%     ylabel('emissions','FontSize',fontAxis);
%     legendStr = {strPVEmission, strGEEmission, strGridEmission};
%     legend(legendStr,'Location','northeast','FontSize',fontLegend);
%     ylim([0,sum(emission_supply_tou)*2.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'emission_DR_sup_only.eps']);
% % 2. Demand only
% % 2.a. Capacity
%     figure;
%     bar([capacity_demand;capacity_demand_tou; capacity_demand_ibr; ...
%         capacity_demand_cpp;capacity_demand_sr;capacity_demand_BiMarkets],groupBarWidth);
%     ylabel('capacity (kW)','FontSize',fontAxis);    
%     legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
%     legend(legendStr,'Location','northwest','FontSize',fontLegend);
%     ylim([0,max(capacity_demand)*1.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'capacity_Demand_DR.eps']);    
%     
% % 2.b. Cost
%     figure;
%     bar([cost_demand; cost_demand_tou; cost_demand_ibr; cost_demand_cpp; ... 
%         cost_demand_sr; cost_demand_BiMarkets],barWidth,'stacked');
%     ylabel('annual expenditure ($)','FontSize',fontAxis);
%     legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
%     legend(legendStr,'Location','northeast','FontSize',fontLegend);
%     ylim([0,sum(cost_demand)*2.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'cost_Demand_DR.eps']);
%     
% % 2.c. Emission
%     figure;
%     bar([emission_demand;emission_demand_tou;emission_demand_ibr; ...
%         emission_demand_cpp;emission_demand_sr; emission_demand_BiMarkets ],barWidth/2,'stacked');
%     ylabel('emissions (gram)','FontSize',fontAxis);
%     legendStr = {strPVEmission, strGEEmission, strGridEmission};
%     legend(legendStr,'Location','northeast','FontSize',fontLegend);
%     ylim([0,sum(emission_demand_tou)*2.5]);
%     xlim([0.5,6.5]);
%     set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
%     set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
%     print ('-depsc', [FIG_PATH 'emission_Demand_DR.eps']);
 

% 3. Integrated
% 3.a. Capacity
    figure;
    bar([capacity_int;capacity_int_tou; capacity_int_ibr; ...
        capacity_int_cpp;capacity_int_sr;capacity_int_BiMarkets],groupBarWidth);
    ylabel('capacity (kW)','FontSize',fontAxis);    
    legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
    legend(legendStr,'Location','southeast','FontSize',fontLegend);
    ylim([0,max(capacity_int_ibr)*1.05]);
    xlim([0.5,6.5]);
    set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', [FIG_PATH 'capacity_int_DR.eps']);    
    
% 3.b. Cost
    figure;
    bar([cost_int;cost_int_tou;cost_int_ibr;cost_int_cpp;cost_int_sr; ...
        cost_int_BiMarkets]*0,barWidth,'stacked');
    ylabel('annual expenditure ($)','FontSize',fontAxis);
    legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
%     legend(legendStr,'Location','southwest','FontSize',fontLegend);
    ylim([0,sum(cost_int_ibr)*1.05]);
    xlim([0.5,6.5]);
    set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', [FIG_PATH 'cost_int_DR.eps']);
    
% 3.c. Emission
    figure;
    bar([emission_int;emission_int_tou;emission_int_ibr; ...
        emission_int_cpp;emission_int_sr; emission_int_BiMarkets ]/1000,barWidth,'stacked');
    ylabel('emissions (kg)','FontSize',fontAxis);
    legendStr = {strPVEmission, strGEEmission, strGridEmission};
    legend(legendStr,'Location','north','FontSize',fontLegend);
    ylim([0,sum(emission_int_BiMarkets)/1000*1.05]);
    xlim([0.5,6.5]);
    set(gca,'xticklabel',{strNo_DR,str_tou,str_ibr,str_cpp,str_sr,str_biMarkets},'FontSize',fontAxis);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', [FIG_PATH 'emission_int_DR.eps']);

%% Capacity and Operation Analysis
figure_settings;
load('results/houston_payback_seasonal_0213.mat');

figure;
bar([cost_grid;cost_supply;cost_demand;cost_int],barWidth/2,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
ylim([0,sum(cost_grid)*2]);
xlim([0.5,4.5]);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_comp.eps']);

figure;
bar([capacity_grid;capacity_supply;capacity_demand;capacity_int],groupBarWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(capacity_demand)*1.5]);
set(gca,'xticklabel',{gridStr,supplyStr,demandStr,integratedStr},'FontSize',fontAxis);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_comp.eps']);

%% grid prices 
figure_settings;
load('results/TCO_grid_electricity_0115.mat');

figure;
xArray = GP_Array;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,GP_Array,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel('electricity price ($/kWh)','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northeast','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.1]);
xlim([0,max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_electricity.eps']);

figure;
plot(GP_Array,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel('electricity price ($/kWh)','FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([0,max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_electricity.eps']);

%% gas prices
figure_settings;
load('results/houston_payback_gas.mat');

figure;
bar(gasPrices,cost,barWidth,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
xlabel('gas price ($/kWh)','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(sum(cost,2))*1.1]);
% xlim([0,0.11]);
xlim([0.01,0.13]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print('-depsc', [FIG_PATH 'cost_gas.eps']);

figure;
xArray = gasPrices;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,xArray,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel('gas price ($/kWh)','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.1]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_gas.eps']);

figure;
plot(xArray,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel('gas price  ($/kWh)','FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([0,max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_gas.eps']);

%% PV Capacity Factor
figure_settings;
load('results/houston_payback_pv.mat');

figure;
bar(CF,cost,barWidth,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
xlabel('PV capacity factor','FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','west','FontSize',fontLegend);
ylim([0,max(sum(cost,2))*1.1]);
xlim([0.11,0.42]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_PV.eps']);

figure;
xArray = CF;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,xArray,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel('PV capacity factor','FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.7]);
xlim([0.13,0.4]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_PV.eps']);

figure;
plot(xArray,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel('PV capacity factor','FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([min(xArray),max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_PV.eps']);


%% Ratio of flexible workload
figure_settings;
load('results/houston_payback_ratio.mat');


xArray = ratio;
xlabelStr = 'Flexible workload ratio';

figure;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,xArray,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','northwest','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.05]);
% xlim([0.13,0.4]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_batch.eps']);

figure;
bar(xArray,cost,barWidth,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
ylim([0,max(sum(cost,2))*1.05]);
xlim([-0.1 1.1]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_batch.eps']);

figure;
plot(xArray,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel('capacity (kW)','FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([min(xArray),max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_batch.eps']);

%% Non-flexible workload
figure_settings;
load('results/houston_payback_PMR.mat');
xArray = PMR;
xlabelStr = 'Peak-to-Mean Ratio';

figure;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,xArray,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','west','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.1]);
xlim([1.3,2.8]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_PMR.eps']);
% saveas(gcf,'results\capacity_comp_pb_gp.fig');

figure;
bar(xArray,cost,barWidth,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','north','FontSize',fontLegend);
ylim([0,max(sum(cost,2))*1.1]);
xlim([1.2,2.9]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_PMR.eps']);
% saveas(gcf,'results\expenditure_comp_pb_gp.fig');

figure;
plot(xArray,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel(xlabelStr,'FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([min(xArray),max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_PMR.eps']);


%% Consolidation
figure_settings;
load('results/houston_payback_con.mat');

xArray = CON;
xlabelStr = 'maximum utlization after consolidation';

figure;
plot(xArray,capacity(:,1),PVLineStyle,xArray,capacity(:,2),GELineStype,xArray,capacity(:,3),gridLineStype, 'linewidth',lineWidth);
ylabel('capacity (kW)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVCapacityStr, GECapacityStr, gridCapacityStr};
legend(legendStr,'Location','west','FontSize',fontLegend);
ylim([0,max(max(capacity))*1.05]);
xlim([min(CON),max(CON)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'capacity_con.eps']);
% saveas(gcf,'results\capacity_comp_pb_gp.fig');

figure;
bar(xArray,cost,barWidth,'stacked');
ylabel('annual expenditure ($)','FontSize',fontAxis);
xlabel(xlabelStr,'FontSize',fontAxis);
legendStr = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill};
legend(legendStr,'Location','southeast','FontSize',fontLegend);
ylim([0,max(sum(cost,2))*1.01]);
xlim([min(CON)-0.05,max(CON)+0.05]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'cost_con.eps']);
% saveas(gcf,'results\expenditure_comp_pb_gp.fig');

figure;
plot(xArray,payback,paybackLineStyle, 'linewidth',lineWidth)
xlabel(xlabelStr,'FontSize',fontAxis);
ylabel('pay back period (years)','FontSize',fontAxis);
ylim([0,max_payback+1]);
xlim([min(xArray),max(xArray)]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
print ('-depsc', [FIG_PATH 'payback_con.eps']);

%%
close all;