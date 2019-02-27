
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the common settings.
init_common_settings_allox;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
nz = 0;
pd = 'normal';

alpha = 10;
beta  = 10;
ITC_gpu = ITC*5;
ITR_gpu = ITR;

% BS_Array = 0*BS_Array;
ds = 0; % scheduling batch jobs.
% RATIOS = [1 2 4 8 16 32];
% RATIOS = [1 2 4 10] ;
RATIOS = [4] ;
a_Array_gpu_org = a_Array_gpu;
progressbar(0);
for iter=1:length(RATIOS)
%     addPredictionErrors;
    CPU_GPU_RATIO = RATIOS(iter);
    a_Array_gpu = a_Array_gpu_org/ mean(a_gpu) * (mean(a) / CPU_GPU_RATIO);
    mean(a_Array_gpu)
    BS_Array_gpu = BS_Array_gpu / mean(BS_Array_gpu) * (mean(BS_Array)) / CPU_GPU_RATIO;
    mean(BS_Array_gpu)
    %% peak provision: interchange = 0
    interchange=0;
    [cost_trad, capacity_trad, emission_trad] = houston_grid_longterm_gpu(interchange, 1,nz,ds,1, ...
        T, N_y, ...
        IP, PP, OP, ...
        ITR,ITC, ...
        ITR_gpu,ITC_gpu,alpha,beta, ... % gpu
        CO2_grid, a_Array, a_Array_gpu ...
        ,au, ...
        con, BS_Array, BS_Array_gpu ...
        ,A_Array, A_Array_gpu, ...
        S,E,bu,PUE_low,R(1,:), ...
        RR,RC_Array ...
        ,RE,SR, ...
        CRR, ...
        CRC_Array ...
        ,CRE,P, GP_Array ...
        ,RP,BP,plot_option,[fig_path 'nz-int_year']);
    costTrads{iter} = cost_trad;  
    capacityTrads{iter} = capacity_trad; 
    emissionTrads{iter} = emission_trad; 
    
    progressbar((iter*2-1)/(length(RATIOS)*2));
    %% integrated with interchange = 1
    interchange = 1;
    [cost_int, capacity_int, emission_int] = houston_grid_longterm_gpu(interchange, 1,nz,ds,1, ...
        T, N_y, ...
        IP, PP, OP, ...
        ITR,ITC, ...
        ITR_gpu,ITC_gpu,alpha,beta, ... % gpu
        CO2_grid, a_Array, a_Array_gpu ...
        ,au, ...
        con, BS_Array, BS_Array_gpu ...
        ,A_Array, A_Array_gpu, ...
        S,E,bu,PUE_low,R(1,:), ...
        RR,RC_Array ...
        ,RE,SR, ...
        CRR, ...
        CRC_Array ...
        ,CRE,P, GP_Array ...
        ,RP,BP,plot_option,[fig_path 'nz-int_year']);
    % realize the installed capacity & optimize the demand
    costInts{iter} = cost_int;  
    capacityInts{iter} = capacity_int; 
    emissionInts{iter} = emission_int; 
    
    progressbar(iter*2/(length(RATIOS)*2));
end

%%
if OFFICIAL
  fileName = strcat('results/TCO_comparison_err_gpu_7years','.mat');
else
  fileName = strcat('results/TCO_comparison_err_gpu_1year','.mat');
end
% fileName = strcat('results/TCO_comparison_err_gpu','.mat');
save(fileName);