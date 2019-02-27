%% long-term TCO_comparison
common_settings;
load('results/TCO_comparison_err_gpu_7years.mat'); expId=5;
% load('results/TCO_comparison_err_gpu_32.mat'); expId=1;
% load('results/TCO_comparison_err_gpu.mat'); expId=1; % for 16
plots = [1 1];
figIdx = 0;
is_printed = 1;
figure_size = figSizeTwothirdCol;
strTraditional = 'traditional';
strOptimal = 'proposed';
%% 
if plots(1)
  figure;
  for expId = 1:length(RATIOS)
    cost_trad =  sum(costTrads{expId}');
    cost_int =  sum(costInts{expId}');  
    costSavings(expId) = (cost_trad(7) - cost_int(7)) / cost_trad(7) *100;
  end 
%   RATIOS = [1 4 16 32];
%   costSavings = [1 2 4.6 15];
%   RATIOS = [1 2 4 10];
%   costSavings = [1 2 3 4.6];

  plot(RATIOS, costSavings, 'LineWidth', LineWidth);
  ylabel('GPU cost saving (%)','FontSize', fontAxis);
%   legendStr = {'CPU', 'GPU'};
%   legendStr = {strTraditional,strOptimal};
%   legend(legendStr,'Location','northwest','FontSize',fontLegend);  
  ylim([0,30]);
  xlim([1 max(RATIOS)]);
%   xLabels = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill, 'CPU', 'GPU'};
  xlabel('CPU to GPU ratio','FontSize',fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figure_size, 'PaperUnits', 'inches', 'PaperPosition', figure_size);
  
  if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'prov_savings';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
   end 
  
end

if plots(2)
  figure;
  expId = 4;
  cost_trad =  sum(costTrads{expId}');
  cost_int =  sum(costInts{expId}');
  
  costsTrad = [cost_trad(1)+cost_trad(2) cost_trad(3)+cost_trad(4) cost_trad(5) cost_trad(6) cost_trad(7)];
  costsInt = [cost_int(1)+cost_int(2) cost_int(3)+cost_int(4) cost_int(5) cost_int(6) cost_int(7)];
  bar([costsTrad;costsInt]','grouped');
  ylabel('expenditure ($)','FontSize',fontAxis);
%   legendStr = {'CPU', 'GPU'};
  legendStr = {strTraditional,strOptimal};
  legend(legendStr,'Location','northwest','FontSize',fontLegend);  
%   ylim([0,max(capacity_sum_int)*1.1]);
%   xlim([0.5 2.5]);
%   xLabels = {PVInfraCost, PVOMCost, GEInfraCost, GEOMCost, electricityBill, 'CPU', 'GPU'};
  xLabels = {'PV', 'GE', 'Util', 'CPU', 'GPU'};
  set(gca,'xticklabel',xLabels,'FontSize', fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figure_size, 'PaperUnits', 'inches', 'PaperPosition', figure_size);
  
  if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'prov_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
   end 
  
end
%%

%%
return;
%% convert to PDFs
fig_path = 'figs/';
for i=1:length(fileNames)
    fileName = fileNames{i}
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end