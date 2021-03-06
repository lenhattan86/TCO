
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the common settings.
init_common_settings_3;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
nz = 0;
pd = 'normal';

predErrStd=0.0;
% N_samples = 5;
if predErrStd==0 
    N_samples = 1;
end
progressbar;
for iter=1:N_samples
    addPredictionErrors;

    %% integrated
    [cost_int, capacity_int, emission_int] = houston_grid_longterm_IT(1,nz,1,1, ...
        T, N_y, ...
        IP, PP, OP, ...
        ITR,ITC, ...
        CO2_grid, a_Array ...
        ,au, ...
        con, BS_Array ...
        ,A_Array, ...
        S,E,bu,PUE_low,R(1,:), ...
        RR,RC_Array ...
        ,RE,SR, ...
        CRR, ...
        CRC_Array ...
        ,CRE,P, GP_Array ...
        ,RP,BP,plot_option,[fig_path 'nz-int_year']);
%     % realize the installed capacity & optimize the demand
%     [cost_int_actual,emission_int_actual] = realize_houston_grid_longterm(1,nz,1,1, ...
%         T, N_y, ...
%         IP,PP,OP,IT,CO2_grid,a_Array_real ...
%         ,au, ...
%         con, BS_Array_real ...
%         , A_Array, ...
%         S,E,bu,PUE_low,R_real(1,:), ...
%         RR,RC_Array ...
%         ,RE,SR, ...
%         CRR, ...
%         CRC_Array_real ...
%         ,CRE,P, GP_Array_real ...
%         ,RP,BP,capacity_int);
    
    progressbar(iter/N_samples);
end
%%
fileName = strcat('results/TCO_comparison_err_IT',num2str(predErrStd),'.mat');
save(fileName);