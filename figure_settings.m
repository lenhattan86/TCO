clear all; clc; close all;
FIG_PATH = 'figs/';
% FIG_PATH = '../manuscript/figs/';
%FIG_PATH = 'results/';

gridOnlyStr = 'Grid Only';
supplyOptStr = 'Supply-side Optimization';
demandOptStr = 'Demand-side Optimization';
integratedOptStr = 'Integrated Optimization';

gridStr = 'GRID';
supplyStr = 'SUP';
demandStr = 'DEM';
integratedStr = 'PROP';

gridCapacityStr = 'Grid';
PVCapacityStr = 'PV';
GECapacityStr = 'GE';
ITCapacityStr= 'IT';

PVInfraCost =       'PV-I.';
PVOMCost =          'PV-O';
GEInfraCost =       'GE-I';
GEOMCost =          'GE-O';
electricityBill =   'Util';

strPVEmission =   'PV';
strGEEmission =   'GE';
strGridEmission = 'Grid';

strNo_DR = 'No-DR';
str_tou  = 'ToU';
str_ibr  = 'IBR';
str_cpp  = 'CPP';
str_sr   =  'SR';
str_biMarkets = 'WS';


PVLineStyle = '-r';
GELineStype = '>--';
gridLineStype = 'x-.b';
ITLineStyle = 'o-k';
paybackLineStyle = 'x-b';

gridLine = '-r';
supLine = '>--';
demLine = 'x-.b';
propLine = 'o-k';

worstOfflineLineStyle = '-.';
bestOfflineLineStyle = ':';

lineWidth = 1.5;
barWidth = 0.4;
groupBarWidth = 0.9;

fontAxis = 12;
fontTitle = 12;
fontLegend = 12;
FontSize = 12;

defaultColorSet = get(gcf,'colormap');

colorset = [0   0   1;...
1   0   0;...
0   1   0;...
0   0   0.172413793103448;...
1   0.103448275862069   0.724137931034483;...
1   0.827586206896552   0;...
0   0.344827586206897   0;...
0.517241379310345   0.517241379310345   1;...
0.620689655172414   0.310344827586207   0.275862068965517];

IS_CHANGE_COLOR = 0;

figure_size = [0.0 0 4.0 3.5];
close all;