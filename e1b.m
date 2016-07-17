clc;
clear;

% constant input
D = 7; % total days
TS = 24; % number of time slots in one day


% PV solar traces input
r_raw = load('traces/solar-one-week.csv'); % renewable trace, in kW
NR = 12;% number of points per hour % number of points per hour
CR = 1; % column of the renewable data

% grid electricity price, from real traces
p_raw = load('traces/pa-electric-price.txt'); % electricity price trace
PS = 0.1; % scale to cents/kWh

% hourly temperature, from real traces
t_raw = load('traces/hourly-temperature.txt');
%COP = 6*ones(TS,1);

RAT = 25;
MaxBlowerCFM=10000;
MaxBlowerPower=11;
NB=1;
beta=3;

% capacity cap, from real traces and model
CC = 1.5; % define the capacity as CC*(peak renewable supply)

% server power cost setting
P_idle = 100; % idle power is 100W
P_peak = 200; % peak power is 200W
nCPU = 12; % number of CPU is 12

% interactive workload traces input
a_raw = load('traces/SAPnew/sapTrace.tab');  % interactive workload trace
ap = load('traces/SAP_prediction.txt')';
NI = 12; % number of points per hour
CI = 4; % column of the interactive data
u = 0.2; % highest utilization level given by SLA
NS = 1e2; % scale factor, changed into CPU share
RI = 2;  % total interative workload/total PV supply

% batch workload
N = 2; % total number of batch workload
B = 2000000*ones(1,2); % batch workload demand in CPU shares
R = [0.000042, 0.000012]; % batch workload revenue per 100 CPU share per hour
maxCPU = [5000, 5000]; % maximum CPU numbers of each batch jobs
S = [1,1,1]; % start time
E = [24,24,24]; % end time

% energy storage
ES0 = 0; % initial energy storage
EC = 0; % energy capacity
Loss = 0.95; % convertion rate

% cooling constant
CoC = 5;

% misc
NZ = 1.0; % net-zero ratio
BA = 0.8; % bandwidth constraint

% figure print control
printRenewable = 0; % print renewable supply and interactive workload
printElectricity = 0; % print electricity price
day = 1; % choose the day to print
has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_clear;


% choice: (1-4) means revenue problem
% 1 means optimal solution with dynamic pricing
% 2 means optimal solution with flat price
% 3 means optimal solution with net-zero constraint
% 4 means renewable-supply-oblivious solution 

%choice1 = menu('Pick one','Revenue from batch workload','Finish all batch workload');

%choice = menu('Pick one setting','Perfect predictions','Full integration','Supply integration','No integration','All','PV prediction');

choice = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        INPUT PROCESS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = TS*D; % total time slots of interests
% renewable traces processing
for t = 1:1:T % sum over each hour, in kW
    r(t) = max(0,sum(r_raw(t*NR-NR+1:t*NR)));
end

% interactive workload traces processing
for t = 1:1:T % sum to each hour
    a(t) = sum(a_raw(t*NI-NI+1+NI*TS*22:t*NI+NI*TS*22,CI));
end
a = a*NS; % scale to CPU share
a_n = a/100; % number of CPU
a_r = a/u*(1-u); % available share for batch workload
a = a/u/100/nCPU*(P_idle + (P_peak-P_idle)*u)/1000; % changed into kW
a(1:TS) = a(1:TS)/mean(a(1:TS))*mean(r)*RI;
% predicted interactive workload traces processing
ap = ap*NS; % scale to CPU share
ap_n = ap/100; % number of CPU
ap_r = ap/u*(1-u); % available share for batch workload
ap = ap/u/100/nCPU*(P_idle + (P_peak-P_idle)*u)/1000; % changed into kW
ap = ap/mean(ap)*mean(r)*RI;
%{
figure;
plot(1:TS, a(1:TS),1:TS,ap(1:TS))
xlim([1,24]);
legend('actural workload','predicted worload');
xlabel('time (hour)');
ylabel('SAP workday workload');
legend('predicted workload','actual workload')
title('SAP workday workload prediction')
%}


% electricity price
p = p_raw*PS; % scale to cents/kWh
p = mean(p)*ones(1,T);


% hourly temperature
OAT = (t_raw - 32)/1.8;
y = 4*400;
COP= 2.617 + 2.938*y./((-19.7*(OAT.^(1/2)))*(-1.105e-005*y*y+0.01681*y-14.54)-4);
AirCp = 1.006; % air capacity, KW/degC
MaxBlowerFlowRate = MaxBlowerCFM * 1.16/2119; % in kg/sec

% capacity cap, from real traces and model % grid electricity price, from real traces
%C = max(r(1:T))*CC; 
C = 100;
% batch workload
for d = 1:1:D
    for n = 1:1:N
        A((d-1)*TS+S(n):(d-1)*TS+E(n),n) = ones(E(n)-S(n)+1,1); % availability matrix
    end
end

p_original = p;
if min(p) <= 0
    p = p - min(p);
    R = R - min(p);
end
b_power = zeros(T,N);
b_record = zeros(T,N);
cooling = zeros(T,1);
b_prediction_power = zeros(T,N);
b_prediction_record = zeros(T,N);
cooling_prediction = zeros(T,1);
b_PV_power = zeros(T,N);
b_PV_record = zeros(T,N);
cooling_PV = zeros(T,1);
b_constant_power = zeros(T,N);
b_constant_record = zeros(T,N);
cooling_constant = zeros(T,1);
b_flat_power = zeros(T,N);
b_flat_record = zeros(T,N);
cooling_flat = zeros(T,1);
rp =  load('traces/papvpred-8-25-2011.txt')';

%{
figure;
plot(1:TS, r(1:TS),1:TS,rp(1:TS))
xlim([1,24]);
legend('actural PV generation','predicted PV generation');
xlabel('time (hour)');
ylabel('PV generation (KW)');
title('PV generation prediction')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a(1:TS) = a(1:TS)/mean(a(1:TS))*18.5;
B = ones(1,2)*sum(a(1:TS))*10000; % batch workload demand in CPU shares
for i = 1:1:10
    a(1:TS) = a(1:TS)/mean(a(1:TS))*i*4; 
    for d = 1:1:1
        cvx_begin
            variables b(TS,N) s(TS) n1(TS,N) n2(TS,N) % b is per CPU share
            minimize( p(d*TS-(TS-1):d*TS)*max(0,a(d*TS-(TS-1):d*TS)' + sum(n1')' * (P_peak-P_idle)/nCPU*(1-u)/1000 + sum(n2')' * P_peak/nCPU/1000+(...
                max(0,a(d*TS-(TS-1):d*TS)' + sum(n1')' * (P_peak-P_idle)/nCPU*(1-u)/1000 + sum(n2')' * P_peak/nCPU/1000 - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)...
                ) - r(d*TS-(TS-1):d*TS)') - R*sum(A(d*TS-(TS-1):d*TS,:).*b)' )
            subject to
                b == n1 * (1-u) * 100 + n2 * 100;
                n1 + n2 <= ones(TS,1) * maxCPU;
                sum(n1')' <= a_n(d*TS-(TS-1):d*TS)';
                n1 >= 0;
                n2 >= 0;
                s >= 0;
                s <= 1;
                %c >=  sum(b')'*(P_peak-P_idle)/nCPU/(100*(1-u))/1000 + max(0,sum(b')-a_r(d*TS-(TS-1):d*TS))'/nCPU*(P_peak/100-(P_peak-P_idle)/(100*(1-u)))/1000 + a(d*TS-(TS-1):d*TS)';
                sum(b)' <= B';
                %c-r(d*TS-(TS-1):d*TS)' <= C*BA*ones(TS,1); % bandwidth constraint
                a(d*TS-(TS-1):d*TS)' + sum(n1')' * (P_peak-P_idle)/nCPU*(1-u)/1000 + sum(n2')' * P_peak/nCPU/1000 <= C;
                %c >= min(C*ones(TS,1), r(d*TS-(TS-1):d*TS)');
                b >= 0;
        cvx_end
        if strcmp(cvx_status, 'Infeasible') == 1
            'Infeasible'
        end
        a_record(1:T) = a;
        n1_record(d*TS-(TS-1):d*TS,:) = n1;
        n2_record(d*TS-(TS-1):d*TS,:) = n2;
        b_record(d*TS-(TS-1):d*TS,:) = b;
        b_power(d*TS-(TS-1):d*TS,:) = n1 * (P_peak-P_idle)/nCPU*(1-u)/1000 + n2 * P_peak/nCPU/1000;
        cooling(d*TS-(TS-1):d*TS) = max(0,a(d*TS-(TS-1):d*TS)' + sum(n1')' * (P_peak-P_idle)/nCPU*(1-u)/1000 + sum(n2')' * P_peak/nCPU/1000 - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
        c_record(d*TS-(TS-1):d*TS) = sum(b_power(d*TS-(TS-1):d*TS,:)')+ a(d*TS-(TS-1):d*TS);
        value_record(d) = cvx_optval;
        consumption(d) = sum(max(0,c_record(d*TS-(TS-1):d*TS)+cooling(d*TS-(TS-1):d*TS)'-r(d*TS-(TS-1):d*TS)));
        renewable(d) = sum(r(d*TS-(TS-1):d*TS)-max(0,r(d*TS-(TS-1):d*TS)-c_record(d*TS-(TS-1):d*TS)-cooling(d*TS-(TS-1):d*TS)'));
        sum(sum(b_power(1:TS,:)))/sum(a(1:TS));
    end
    eff(i) = 5000*(consumption)/sum(sum(b_record));
    ratio(i) = sum(sum(b_power))/(sum(sum(b_power))+sum(a(1:TS))) 
    
    
    b_flat_record = zeros(T,N);
    b_flat_power = zeros(T,N);
    for j = 1:1:N
        b_flat_record((d-1)*TS+1:(d-1)*TS+TS,j) = sum(b_record(1:24,j))/TS;
        b_flat_power((d-1)*TS+1:(d-1)*TS+TS,j) = sum(b_power(1:24,j))/TS;
    end
    c_tmp_record(d*TS-(TS-1):d*TS) = a(d*TS-(TS-1):d*TS)'+sum(b_flat_power(d*TS-(TS-1):d*TS,:)')';
    for d = 1:1:1
        cvx_begin
            variables s(TS) % b is per CPU share
            minimize( p(d*TS-(TS-1):d*TS)*max(0,c_tmp_record(d*TS-(TS-1):d*TS)'+(max(0,c_tmp_record(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)) - r(d*TS-(TS-1):d*TS)'))
            subject to
                s >= 0;
                s <= 1;
        cvx_end
        c_flat_record(d*TS-(TS-1):d*TS) = c_tmp_record(d*TS-(TS-1):d*TS)';
        cooling_flat(d*TS-(TS-1):d*TS) = (max(0,c_tmp_record(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta))';
        consumption_flat(d) = sum(max(0,c_flat_record(d*TS-(TS-1):d*TS)+cooling_flat(d*TS-(TS-1):d*TS)'-r(d*TS-(TS-1):d*TS)));
        renewable_flat(d) = sum(r(d*TS-(TS-1):d*TS)-max(0,r(d*TS-(TS-1):d*TS)-c_flat_record(d*TS-(TS-1):d*TS)-cooling_flat(d*TS-(TS-1):d*TS)'));
    end
    eff_flat(i) = 5000*(consumption_flat)/sum(sum(b_flat_record));
    ratio_flat(i) = sum(sum(b_flat_power))/(sum(sum(b_flat_power))+sum(a(1:TS))) 
    
    
end
figure;
plot(ratio_flat,1-eff./eff_flat,'b')
%ylim([0,1.2*max(compare2)]);
%title('Efficiency comparisons');
xlabel('percentage of batch jobs')
ylabel('energy efficiency improvement');
xlim([0.2,0.8]);
%legend('Flat','Optimal')
%set(gca,'xticklabel',{'50%', '30%','30% Flat'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
print ('-depsc', 'results/e1-comp2.eps');

    
if choice == 5
    compare1 = [renewable,renewable_prediction,renewable_flat;consumption,consumption_prediction,consumption_flat];
    figure;
    bar(compare1',0.3,'stacked')
    ylabel('power (KWh)');
    ylim([0,1500]);
    title('Power consumption comparisons');
    legend('renewable','non-renewable');
    set(gca,'xticklabel',{'50%', '30%','30%+flat'});
    compare2 = 1000*[(consumption)/sum(sum(b_record)),(consumption_prediction)/sum(sum(b_prediction_record)),(consumption_flat)/sum(sum(b_flat_record))];
    %compare2 = 1./compare2;
    figure;
    bar(compare2',0.3)
    %ylim([0,1.2*max(compare2)]);
    title('Efficiency comparisons');
    ylabel('kWh grid power/jobs');
    %legend('grid power','renewable')
    set(gca,'xticklabel',{'50%', '30%','30% Flat'});
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.2]);
    print ('-depsc', 'results/e1-comp2.eps');
    saveas(gcf,'results/e1-comp2.fig')
end
% time | interative demand in kW | interative CPU | batch workload demand
% in kW | batach workload shared CPU | batch workload own CPU | total power
% | renewable supply in kW 


%csvwrite('results\flat.csv',[(1:TS);a(1:TS);a_n(1:TS);b_power(1:TS,:)'; n1(1:TS,:)'; n2(1:TS,:)'; c_record(1:TS); r(1:TS)]',1,0);
    
