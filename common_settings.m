clear; close all; clc;

fontSize=10;
fontAxis = fontSize;
fontTitle = fontSize;
fontLegend = fontSize;
LineWidth = 1.5;
FontSize = fontSize;
is_printed = true;
axisWidth = 1.5;
%%

legendSize = [0.0 0 5 0.4];
figSizeOneCol = [0.0 0 5 3];
figSizeOneColHaflRow = [1 1 1 0.5].* figSizeOneCol;
figSizeTwoCol = 2*figSizeOneCol;
figSizeOneThirdCol = 1/3*figSizeOneCol;
figSizeHalfCol = 1/2*figSizeOneCol;
figSizeTwothirdCol = 2/3*figSizeOneCol;
figSizeFourFifthCol = 4/5*figSizeOneCol;

%%

barLineWidth=0;
groupBarSize = 0.9;
barSize = 0.5;

%%

figIdx=0;

LOCAL_FIG = 'figs/';

%PS_CMD_FORMAT='ps2pdf -dEmbedAllFonts#true -dSubsetFonts#true -dEPSCrop#true -dPDFSETTINGS#/prepress %s %s';
PS_CMD_FORMAT='ps2pdf -dEmbedAllFonts#true -dSubsetFonts#true -dEPSCrop#false -dPDFSETTINGS#/prepress %s %s';

% fig_path = ['figs/'];
fig_path = 'figs/';

%%

strUser1 = 'User 1';
strUser2 = 'User 2';
strUser3 = 'User 3';
strUser4 = 'User 4';
strUnalloc = 'unallocated';

strES = 'ES';
strDRF = 'DRF';
strFDRF = 'FDRF';
strMP = 'MP';
strMSR = 'MSR';
strPricing = 'Pricing';

strGPU = 'GPU';
strCPU = 'CPU';
strMemory= 'mem.';

strJobCompleted = 'completed jobs';

strNormCapacity='norm. capacity';


strEstimationErr = 'std. of estimation errors [%]';
strFactorImprove = 'factor of improvement';
strAvgComplTime = 'avg. compl. (secs)';

strMethods='methods';

%% line specs
lineProposed = '-';
lineStrict = '+:';
lineDRF = '-.';
lineDRFW = '--';

lineBB = '-';
lineTPCDS = '--';
lineTPCH = '-.';
workloadLineStyles = {lineBB, lineTPCDS, lineTPCH};

%%

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