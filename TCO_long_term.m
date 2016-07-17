% based on TCO_grid_electricity_0115.m
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
load_exist = 0;
    
%% Load the common settings.
init_common_settings_3;
genAnnualTrend;

%% additional parameters by Tan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
discount_factor = 0;
demand_opt = 1;
nz = 0;

[cost,capacity, emissions, total_demand, total_supply] = houston_grid_longterm(1,nz,1,1, ...
    T, N_y, ...
    IP,PP,OP,IT,CO2_grid,a_Array ...
    ,au, ...
    con, BS_Array ...
    ,A_Array, ...
    S,E,bu,PUE_low,R(1,:), ...
    RR,RC_Array ...
    ,RE,SR, ...
    CRR, ...
    CRC_Array ...
    ,CRE,P, GP_Array ...
    ,RP,BP,plot_option,'figs/temp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
save('results/TCO_long_term.mat');