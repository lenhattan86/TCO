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
TS = 24; % number of slots per day
ND = 12; % number of days
T = TS*ND; % planning period length, in hours
IP = 0.04; % idle power consumption of a server, in kW
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
for c = 1:1:5 
    DC_power = 200+c*100; % Data center power capacity, in kW
    OP = 0.8; % over-provisioning ratio of IT and power
    IT = DC_power*OP;% Data center IT capacity, in kW
    con = 0.9; % maximum utilization after consolidation
    plot_PV = 0; % plot PV traces
    plot_demand = 0; % plot power demand

    % interactive workloads
    IN = 1; % number of interactive workloads
    IL = {'traces\SAPnew\sapTrace2.tab'};
    a = zeros(IN,T);
    au = 0.4; % average utilization of interactive workloads
    IPU = 0.7; % interactive workload peak utilization
    IM = IT*au*IPU; % peak of interactive workload consumption, in servers
    PMR = 2.5; % peak to mean ratio, should be smaller than current PMR
    for i = 1:1:IN
        a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
    end


    % batch job
    BN = 1; % average number of batch job arrivals per timeslot
    BM = 0.25; % batch job ratio, compared with interactive workload
    [A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties
    bu = 1; % average utilization of batch job when running alone.

    % uncontrolled renewable generation, e.g., PV, wind
    RN = 1; % number of renewable generation
    RL = {'traces\solar-one-week-07012012.csv'};
    RCF = load('traces\PV_CF.csv');
    R = zeros(IN,T);
    RP = 0.3; % capacity factor
    RR = [0, 500]; % renewable capacity range, in kW
    for i = 1:1:RN
        R_raw = trace_process(char(RL(i)), TS*7, 12, 1, 4, TS);
        R(i,:) = reshape(R_raw'/mean(R_raw) * RCF,1,T)/100;
    end
    RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
    RE = 32; % emmision, in gram CO2 eq

    % plot PV traces
    if plot_PV
        figure;
        plot(1/TS:1/TS:ND, R(1,:)*100,'r-','LineWidth',2)
        ylabel('generation (kW)')
        ylim([0,100])
        set(gca,'XTick',[0.5:1:11.5], 'FontSize', 8);
        set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
    end

    % controlled renewable generation, e.g., gas engine
    CRN = 1; % number of controlled renewable generation
    SR = 1/6; % storage requirement for 1kW, in kWh
    CRR = [0,1400]; % controled renewable capacity range, in kW
    NG_price = load('traces\NG_price.csv');
    %CRC = [1000/20+SR*200/4, 0.005+reshape(ones(TS,1) * NG_price',1,T)]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
    CRC = [1000/20+SR*200/4, 0.005+0.10*ones(1,T)];
    CRE = 443; % emmision, in gram CO2 eq

    % cooling efficiency
    PUE = 1.2*ones(1,T);

    % price setting [power, cooling, PV, IT]
    P = [600, 800, 3200/20, 1000];
    GP = reshape(ones(TS,1)*load('traces\electricity_price.csv')',1,T); % grid electricity price
    %GP = 0.10*ones(1,T);
    BP = 0.04*ones(1,T); % sell back price
    RP = 0.14*ones(1,T); % renewable electricity price

    % energy storage loss
    Loss = [0.99, 1];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          SOLVER                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_option = 0;
    discount_factor = 0;
    max_payback = 20;
    payback_plot = 0;

    % grid only
    [cost_grid,capacity_grid,emission_grid] = houston_grid(0,0,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[0,0],RC,RE,SR,[0,0],CRC,CRE,P,GP,RP,BP,max(plot_option + plot_demand*2),'results\tmp');
    yearly_cost_grid = sum(cost_grid([2,4,5]));
    total_cost_grid = yearly_cost_grid*max_payback

    % net-zero, supply given
    plot_option = 0;
    payback_plot = 0;
    for i = 1:1:10
        PV_capacity(i) = i*100;
        [cost_PV(c,i,:),capacity_PV(c,i,:),emission_PV(c,i,:),ratio_PV(c,i,:)] = houston_grid(0,1,0,0,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV_capacity(i),PV_capacity(i)],RC,RE,SR,[250,250],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
        yearly_cost_PV(c) = sum(cost_PV(c,i,[2,4,5]));
        infra_cost_PV(c) = sum(cost_PV(c,i,[1,3]))*max_payback;
        total_cost_PV(c,i) = yearly_cost_PV(c)*max_payback + infra_cost_PV(c)
        cost_reduction_PV(c,i) = 1- total_cost_PV(c,i)/total_cost_grid
        [payback_PV(c,i),npv_PV(c)] = payback_cal(infra_cost_PV(c),(yearly_cost_grid-yearly_cost_PV(c))*ones(1,max_payback),discount_factor,max_payback,payback_plot)

        [cost_PVint(c,i,:),capacity_PVint(c,i,:),emission_PVint(c,i,:),ratio_PVint(c,i,:)] = houston_grid(1,1,1,1,T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,bu,PUE,R(1,:),[PV_capacity(i),PV_capacity(i)],RC,RE,SR,[250,250],CRC,CRE,P,GP,RP,BP,plot_option,'results\tmp');
        yearly_cost_PVint(c) = sum(cost_PVint(c,i,[2,4,5]));
        infra_cost_PVint(c) = sum(cost_PVint(c,i,[1,3]))*max_payback;
        total_cost_PVint(c,i) = yearly_cost_PVint(c)*max_payback + infra_cost_PVint(c)
        cost_reduction_PVint(c,i) = 1- total_cost_PVint(c,i)/total_cost_grid
        [payback_PVint(c,i),npv_PVint(c)] = payback_cal(infra_cost_PVint(c),(yearly_cost_grid-yearly_cost_PVint(c))*ones(1,max_payback),discount_factor,max_payback,payback_plot)
    end
end

figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,total_cost_PV(c,:));
end
box;
ylabel('total expenditure ($)');
xlabel('PV capacity');
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
xlim([100,1000]);
ylim([min(min(total_cost_PV*0.9)),max(max(total_cost_PV*1.1))]);

figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,ratio_PV(c,:,1));
end
box;
xlabel('PV capacity');
ylabel('demand powered by PV');
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
%legend('Total power 300kW','Total power 400kW','Total power 500kW','Total power 600kW','Total power 700kW','Location','Southeast')
xlim([100,1000]);
%ylim([0,1]);


figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,payback_PV(c,:));
end
box
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([100,1000]);
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
ylim([0,20]);


figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,total_cost_PVint(c,:));
end
box;
ylabel('total expenditure ($)');
xlabel('PV capacity');
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
xlim([100,1000]);
ylim([min(min(total_cost_PVint*0.9)),max(max(total_cost_PVint*1.1))]);

figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,ratio_PVint(c,:,1));
end
box;
xlabel('PV capacity');
ylabel('demand powered by PV');
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
xlim([100,1000]);
ylim([0,1]);


figure;
for c = 1:1:5
    hold all
    plot(PV_capacity,payback_PVint(c,:));
end
box
ylabel('payback period (year)');
xlabel('PV capacity');
xlim([100,1000]);
legend('Average IT power 108kW','Average IT power 144kW','Average IT power 180kW','Average IT power 216kW','Average IT power 252kW','Location','Southeast')
ylim([0,20]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save('results\houston_payback_seasonal_IT.mat');
