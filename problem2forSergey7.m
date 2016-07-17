clear;
clc;

choice = menu('Pick one setting','Optimal','Night','All');
switch choice
    case 1,
        scheduling = 1;
    case 2, 
        scheduling = 2;
    case 3,
        scheduling = 0;
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constant input
D = 1; % total days
TS = 24; % number of time slots in one day
T = TS*D; % total time slots of interests
TR = 1000; % mean of transactional workload
PV = 2000; % mean of pv solar supply
N = 24; % total number of batch workload, from model
SN = 6; % number of Servers
CS = 180; % CPU shares per CPU
t_raw = load('traces/hourly-temperature.txt');
NO = 24/TS;
OAT = zeros(1,T);
for t = 1:1:T
    OAT(t) = (mean(t_raw(t*NO-(NO-1):t*NO)) - 32)/1.8;
end
y = 4*400;
COP= 4.246 + 4766*y./(OAT(1:T).^(1/2)*(y*y-1548*y+8.389e5)+5.863e5 );
PUE = 1+1./COP;
%PUE = [1.16,1.17,1.16,1.20,1.22,1.22,1.24,1.26,1.35,1.32,1.25,1.30,1.29,...
%    1.35,1.32,1.40,1.40,1.25,1.29,1.30,1.28,1.29,1.18,1.13];

% transactional workload traces
%L = load('traces\VDR\VDR.txt'); % VDR workload trace for transactional workload
L1 = load('traces/7/nc_kvm_7_web.tab');  % Web workload trace for transactional workload
L2 = load('traces/7/nc_kvm_7_app.tab');  % App workload trace for transactional workload
L3 = load('traces/7/nc_kvm_7_db.tab');  % DB workload trace for transactional workload
LN = 360*24/TS; % 360 means 24hour version, 180 means 12hour version
L = L1(1:LN*TS,6)+L2(1:LN*TS,6)+L3(1:LN*TS,6); % total interactive workload
for t = 1:1:T % average to each hour
    a(t) = mean(L(t*LN-LN+1:t*LN)); % in CPU shares
end
IN = 8; % number of interactive workload instancs
a = a*IN;


% PV solar traces
PAPV = load('traces/solar-one-week.csv'); % PV solar traces
NR = 12*24/TS;
for t = 1:1:T % average to each hour
    r(t) = max(0,sum(PAPV(t*NR-(NR-1):t*NR,1)));
end
Peak = 230;
Idle = 140;
r = r/max(r)*SN*Peak*0.5; % change the peak to 400W;
r_p = 1.4*r; % save the renewable supply in Watt
r = r/Peak*CS*N;
r_r = (r_p-ceil(a/N/CS)*Idle.*PUE + a/N/CS*(Peak-Idle).*PUE)/Peak*CS*N;

% renewable supply in CPU shares, here 2 servers * N CPU * 1.8GHz




%debt = max(0,a_p.*PUE-r);
%r = (r-ceil(r/Peak)*Idle)/(Peak-Idle)*CS*N;

% a = a/mean(a)*mean(r)*0.3; % scale the average to 30% of the renewable supply


%{
figure;
plot(1:T,r,'b',1:T,a,'r');
xlabel('hour');
ylabel('power (CPU shares)');
xlim([1,T]);
legend('renewable supply','transactional workload');
set(gca,'XTick',[0:T/4:T]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', 'results/power.eps');
saveas(gcf,'results/power.fig')
%}

% capacity cap, from real traces and model % grid electricity price, from real traces
C = max(r(1:TS));

% grid electricity price, from real traces
%p = load('traces\pa-electric-price.txt');
p = 0.1*ones(1,TS); % flat electricity price @ 10cents/kWh

%{
figure;
plot(p);
xlabel('hour');
ylabel('electricity price');
xlim([0,T]);
set(gca,'XTick',[0:T/4:T]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', 'results/electricity-price.eps');
saveas(gcf,'results/electricity-price.fig')
%}

p_original = p;
if min(p) <= 0
    p = p - min(p);
end
pb = 0*p; % sell-back price, 0 means no sell-back scheme
pp = 0; % peak demand charging, 0 means no peak demand charging

% batch workload
B = CS*[8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36,8*61/36,8*75/36,1*280*2.4/1.8/36];% batch workload demand
B = B*360/LN;
MP = CS*[8,8,1,8,8,1,8,8,1,8,8,1,8,8,1,8,8,1,8,8,1,8,8,1];% maximum CPU shares

A = zeros(T, N); % availability

S = ones(1,N); % start time
E = TS*ones(1,N); % end time
for d = 1:1:D
    for n = 1:1:N
        A((d-1)*TS+S(n):(d-1)*TS+E(n),n) = ones(E(n)-S(n)+1,1);
    end
end
for i = 1:1:N
    R(i) = 0.4 - 0.01*i;
end
RAT = 35;
OAS = 50; % scale to be suitable for the input
MaxBlowerCFM=10000*OAS;
MaxBlowerPower=11*OAS;
NB=1;
beta=3;
AirCp = 1.006; % air capacity, KW/degC
MaxBlowerFlowRate = MaxBlowerCFM * 1.16/2119; % in kg/sec

% energy storage
ES0 = 0; % initial energy storage
EC = 0; % energy capacity
Loss = 0.95; % loss rate
Conv = 0.9; % convertion rate


% the day of figure
day = 1;

% output file
s1 = 'results/netzero-exp.csv';
s2 = 'results/night-exp.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if scheduling == 1 || scheduling == 0
    b_record = zeros(T,N);
    MF = zeros(TS,1)/Peak*CS*N;
    for d = 1:1:D
        cvx_begin
            variables b(TS,N) c(TS) s(TS) ES(TS) EU(TS)
            minimize( pb(d*TS-(TS-1):d*TS)*(a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)+...
                (p(d*TS-(TS-1):d*TS)-pb(d*TS-(TS-1):d*TS))*max(0,a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)+...
                pp*max((max(0,a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)))- R*sum(A(d*TS-(TS-1):d*TS,:).*b)' )
            subject to
                c  >= a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
                %sum(A(d*TS-(TS-1):d*TS,:).*b)' == B';
                sum(b) <= B;
                b <= ones(TS,1)*MP;
                sum(c+MF) <=  sum(r(d*TS-(TS-1):d*TS));
                %c-r(d*TS-(TS-1):d*TS)' <= C/2*ones(TS,1); % bandwidth constraint
                %c >= min(r(d*TS-(TS-1):d*TS)', C); % free renewable
                c <= C;
                b >= 0;
                s >= 0;
                s <= 1;
                ES + EU/Conv == Loss*[ES0;ES(1:TS-1)];
                ES >= 0;
                ES <= EC;
                ES(TS) == ES0;
        cvx_end


        cooling = max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
        cooling1 = max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'));
        cooling2 = NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT');
        total = [a;b';cooling'];
        figure;
        bar(total(:,TS*(day-1)+1:TS*day)','stacked');
        hold on
        plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
        xlabel('hour');
        ylabel('both workload');
        ylim([0,2*C])
        xlim([1,TS]);
        set(gca,'XTick',[0:TS/4:TS], 'FontSize', 8);

        for i = 1:1:TS
            if abs(c(i)) < 1e-5
                c(i) = 0;
            end
            for n = 1:1:N
                if abs(b(i,n)) < 1e-3
                    b(i,n) = 0;
                end
            end
        end
        c_record(d*TS-(TS-1):d*TS) = a+sum(b');
        b_record(d*TS-(TS-1):d*TS,:) = b;
        optimal(d) = cvx_optval;
        consumption(d*TS-(TS-1):d*TS) = (ceil(c_record(d*TS-(TS-1):d*TS)/CS/N) *Idle + c_record(d*TS-(TS-1):d*TS)/CS/N*(Peak-Idle));
        %consumption(d*TS-(TS-1):d*TS) = c_record(d*TS-(TS-1):d*TS)/CS/N*Peak;
        %cooling(d*TS-(TS-1):d*TS) = consumption(d*TS-(TS-1):d*TS).*(PUE-1);
        brown(d*TS-(TS-1):d*TS) = max(0,consumption(d*TS-(TS-1):d*TS)+cooling'-r_p);
        brown_s = sum(brown);
        %consumption(d*TS-(TS-1):d*TS) = SN *Idle + max(0,c-r(d*TS-(TS-1):d*TS)')/CS/N*(Peak-Idle);
        utilization(d) = 1-sum(max(0,r(d*TS-(TS-1):d*TS)'-c))/sum(r(d*TS-(TS-1):d*TS));
    end
 
    figure;
    a_p = ceil(a/CS/N)*Idle + a/CS/N*(Peak-Idle);
    b_p = consumption - a_p;
    %cooling = cooling/N/CS*Peak;
    cooling = max(0,a_p' + b_p' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak)./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)/CS/N*Peak;
    chiller_power = max(0,a_p' + b_p' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak)./(COP');
    oa_power = NB*MaxBlowerPower * pow_pos(s,beta)/CS/N*Peak;
    oa_capacity = NB*AirCp*MaxBlowerFlowRate * ones(TS,1) .*max(0,RAT' - OAT')/CS/N*Peak;
    chiller = max(0,a_p' + b_p' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak);
    oa = min(consumption',NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak);
    r_p = r_p/sum(r_p)*sum(consumption+cooling');
    total = [a_p;b_p;cooling']/1000; % change to kW
    bar(total','stacked');
    hold on
    plot(1:TS, r_p(TS*(day-1)+1:TS*day)/1000,'-rs', 'LineWidth', 2)
    ylim([0,1.5*max(r_p)/1000]);
    ylabel('power (KW)')
    xlabel('time (0.5 hour)')
    title('Optimal')
    legend('critical workload','non-critical workload','cooling power','renewable supply')
    xlim([1,TS]);
    set(gca,'XTick',[0:TS/4:TS], 'FontSize', 8);
    %{
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    hold on
    bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    xlabel('hour');
    ylabel('cpu shares');
    ylim([0,1.5*C]);
    legend('capacity','renewable supply','interactive','batch')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-opt.eps');
    saveas(gcf,'results/solution-opt.fig')
    %}
    csvwrite(s1,[(1:TS);r_p(1:TS);consumption(1:TS);cooling(1:TS)';chiller_power';oa_power';chiller';oa';oa_capacity';ceil(c/CS/N)';a(1:TS);b_record(1:TS,:)']',1,0);
end
%csvwrite('results/netzero12hour.csv',[(1:TS);r_p(1:TS);consumption(1:TS);a_p;b_p;cooling(1:TS)';ceil(c/CS/N)';a(1:TS);b_record(1:TS,:)']',1,0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          BASELINES                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if scheduling == 2 || scheduling == 0
    b_night_record = zeros(T,N);
    MF = zeros(TS,1)/Peak*CS*N;
    S = ones(1,N); % start time
    E = TS/4*ones(1,N); % end time
    A = zeros(T, N); % availability
    for d = 1:1:D
        for n = 1:1:N
            A((d-1)*TS+S(n):(d-1)*TS+E(n),n) = ones(E(n)-S(n)+1,1);
        end
    end
    for d = 1:1:D
        cvx_begin
            variables b(TS,N) c(TS) s(TS) ES(TS) EU(TS)
            minimize( pb(d*TS-(TS-1):d*TS)*(a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)+...
                (p(d*TS-(TS-1):d*TS)-pb(d*TS-(TS-1):d*TS))*max(0,a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)+...
                pp*max((max(0,a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,max(a(d*TS-(TS-1):d*TS)' + sum(b')',CS*TS) - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)+MF-r(d*TS-(TS-1):d*TS)'-EU)))- R*sum(A(d*TS-(TS-1):d*TS,:).*b)' )
            subject to
                c  >= a(d*TS-(TS-1):d*TS)' + sum(b')'+max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
                %sum(A(d*TS-(TS-1):d*TS,:).*b)' == B';
                sum(b) <= B;
                b <= ones(TS,1)*MP;
                sum(c+MF) <=  sum(r(d*TS-(TS-1):d*TS));
                %c-r(d*TS-(TS-1):d*TS)' <= C/2*ones(TS,1); % bandwidth constraint
                %c >= min(r(d*TS-(TS-1):d*TS)', C); % free renewable
                c <= C;
                b >= 0;
                s >= 0;
                s <= 1;
                ES + EU/Conv == Loss*[ES0;ES(1:TS-1)];
                ES >= 0;
                ES <= EC;
                ES(TS) == ES0;
        cvx_end


        cooling_night = max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
        cooling1_night = max(0,a(d*TS-(TS-1):d*TS)' + sum(b')' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'));
        cooling2_night = NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT');
        total_night = [a;b';cooling_night'];
        figure;
        bar(total_night(:,TS*(day-1)+1:TS*day)','stacked');
        hold on
        plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
        xlabel('hour');
        ylabel('both workload');
        ylim([0,2*C])
        xlim([1,TS]);
        set(gca,'XTick',[0:TS/4:TS], 'FontSize', 8);

        for i = 1:1:TS
            if abs(c(i)) < 1e-5
                c(i) = 0;
            end
            for n = 1:1:N
                if abs(b(i,n)) < 1e-3
                    b(i,n) = 0;
                end
            end
        end
        c_night_record(d*TS-(TS-1):d*TS) = a+sum(b');
        b_night_record(d*TS-(TS-1):d*TS,:) = b;
        optimal_night(d) = cvx_optval;
        consumption_night(d*TS-(TS-1):d*TS) = (ceil(c_night_record(d*TS-(TS-1):d*TS)/CS/N) *Idle + c_night_record(d*TS-(TS-1):d*TS)/CS/N*(Peak-Idle));
        %consumption(d*TS-(TS-1):d*TS) = c_record(d*TS-(TS-1):d*TS)/CS/N*Peak;
        %cooling(d*TS-(TS-1):d*TS) = consumption(d*TS-(TS-1):d*TS).*(PUE-1);
        brown_night(d*TS-(TS-1):d*TS) = max(0,consumption_night(d*TS-(TS-1):d*TS)+cooling_night'-r_p);
        brown_s_night = sum(brown_night);
        %consumption(d*TS-(TS-1):d*TS) = SN *Idle + max(0,c-r(d*TS-(TS-1):d*TS)')/CS/N*(Peak-Idle);
        utilization_night(d) = 1-sum(max(0,r(d*TS-(TS-1):d*TS)'-c))/sum(r(d*TS-(TS-1):d*TS));
    end

    
    figure;
    a_p = ceil(a/CS/N)*Idle + a/CS/N*(Peak-Idle);
    b_p_night = consumption_night - a_p;
    %cooling = cooling/N/CS*Peak;
    cooling_night = max(0,a_p' + b_p_night' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak)./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)/CS/N*Peak;
    chiller_power_night = max(0,a_p' + b_p_night' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak)./(COP');
    oa_power_night = NB*MaxBlowerPower * pow_pos(s,beta)/CS/N*Peak;
    oa_capacity_night = NB*AirCp*MaxBlowerFlowRate * ones(TS,1) .*max(0,RAT' - OAT')/CS/N*Peak;
    chiller_night = max(0,a_p' + b_p_night' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak);
    oa_night = min(consumption_night',NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT')/CS/N*Peak);
    r_p = r_p/sum(r_p)*sum(consumption_night+cooling_night');
    total_night = [a_p;b_p_night;cooling_night']/1000; % change to kW
    bar(total_night','stacked');
    hold on
    plot(1:TS, r_p(TS*(day-1)+1:TS*day)/1000,'-rs', 'LineWidth', 2)
    ylim([0,1.5*max(r_p)/1000]);
    ylabel('power (KW)')
    xlabel('time (0.5 hour)')
    title('Optimal')
    legend('critical workload','non-critical workload','cooling power','renewable supply')
    xlim([1,TS]);
    set(gca,'XTick',[0:TS/4:TS], 'FontSize', 8);
    %{
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    hold on
    bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    xlabel('hour');
    ylabel('cpu shares');
    ylim([0,1.5*C]);
    legend('capacity','renewable supply','interactive','batch')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-opt.eps');
    saveas(gcf,'results/solution-opt.fig')
    %}
    csvwrite(s2,[(1:TS);r_p(1:TS);consumption_night(1:TS);cooling_night(1:TS)';chiller_power_night';oa_power_night';chiller_night';oa_night';oa_capacity_night';ceil(c/CS/N)';a(1:TS);b_night_record(1:TS,:)']',1,0);
end
%csvwrite('results/netzero12hour.csv',[(1:TS);r_p(1:TS);consumption(1:TS);a_p;b_p;cooling(1:TS)';ceil(c/CS/N)';a(1:TS);b_record(1:TS,:)']',1,0);







%{
spare = zeros(1,TS*D);
b_fcfs = zeros(TS*D,N);
b_fcfs_record = zeros(TS*D,N);
c_fcfs_record = zeros(1,TS*D);

if scheduling == 2 || scheduling == 0
for d = 1:1:D
    % baseline 1: FCFS
    spare(d*24-(TS-1):d*24) = 2*24*CS*ones(1,24);   
    for i = 1:1:N
        remain = sum(b_record(:,i));
        current = S(i)+(d-1)*24-1;
        while remain > 0
            current = current +1;
            if current > E(i)+(d-1)*24
                !echo infeasible!
                i
                remain
                %remain = 0;
            end
            if min(remain,MP(i)) <= spare(current)
                b_fcfs(current,i) = min(remain,MP(i));
                spare(current) = spare(current) - min(remain,MP(i));
                remain = max(0,remain - min(remain,MP(i)));
                continue;
            else
                b_fcfs(current,i) = spare(current);
                remain = max(0, remain - spare(current));
                spare(current) = 0;
            end
        end
    end    
    %b_fcfs_record(d*24-(TS-1):d*24,:) = b_fcfs(d*24-(TS-1):d*24,:);
    %fcfs(d) = p(d*24-(TS-1):d*24)*max(0,c_fcfs_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)');
    %consumption_fcfs(d) = sum(max(0,c_fcfs_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)'));
    %utilization_fcfs(d) = 1-sum(max(0,r(d*24-(TS-1):d*24)-c_fcfs_record(d*24-(TS-1):d*24)))/sum(r(d*24-(TS-1):d*24));
end
    c_fcfs(d*24-(TS-1):d*24) = sum(b_fcfs(d*24-(TS-1):d*24,:)');
    consumption_fcfs(d*24-(TS-1):d*24) = ceil(a./PUE/24/CS)*Idle + a./PUE/24/CS*(Peak-Idle)+(ceil(c_fcfs/CS/24) *Idle + c_fcfs/CS/24*(Peak-Idle));
        cvx_begin
            variables s(TS) % b is per CPU share
            minimize( p(d*TS-(TS-1):d*TS)*max(0,consumption_fcfs(d*TS-(TS-1):d*TS)'+(max(0,consumption_fcfs(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)) - r(d*TS-(TS-1):d*TS)'))
            subject to
                s >= 0;
                s <= 1;
        cvx_end
    %consumption(d*24-(TS-1):d*24) = c_record(d*24-(TS-1):d*24)/CS/24*Peak;
    cooling_fcfs(d*24-(TS-1):d*24) = (max(0,consumption_fcfs(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta));
    brown_fcfs(d*24-(TS-1):d*24) = max(0,consumption_fcfs(d*24-(TS-1):d*24)+cooling_fcfs-r_p);
    brown_fcfs_s = sum(brown_fcfs);
    %consumption(d*24-(TS-1):d*24) = SN *Idle + max(0,c-r(d*24-(TS-1):d*24)')/CS/24*(Peak-Idle);
    %utilization(d) = 1-sum(max(0,r(d*24-(TS-1):d*24)'-c))/sum(r(d*24-(TS-1):d*24));

%r_p = r_p/sum(r_p)*sum(consumption+cooling);
end
day = 1;
if scheduling == 2 || scheduling == 0
    
    figure;
    a_fcfs_p = ceil(a./PUE/24/CS)*Idle + a./PUE/24/CS*(Peak-Idle);
    b_fcfs_p = consumption_fcfs - a_p;
    
    total = [a_fcfs_p;b_fcfs_p;cooling_fcfs]/1000;
    bar(total','stacked');
    hold on
    plot(1:TS, r_p(TS*(day-1)+1:TS*day)/1000,'-rs', 'LineWidth', 2)
    ylim([0,1.5 *max(r_p)/1000]);
    ylabel('power (KW)')
    xlabel('time (0.5 hour)')
    title('FCFS');
    legend('critical workload','non-critical workload','cooling power','renewable supply')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
csvwrite('results\netzero12hourfcfs.csv',[(1:TS);r_p(1:TS);consumption_fcfs(1:TS);a_fcfs_p;b_fcfs_p;cooling_fcfs(1:TS);ceil(c_fcfs/CS/24)+1;a(1:TS);b_fcfs(1:TS,:)']',1,0);
    %{
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    hold on
    bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    xlabel('hour');
    ylabel('cpu shares');
    ylim([0,1.5*C]);
    legend('capacity','renewable supply','interactive','batch')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-opt.eps');
    saveas(gcf,'results/solution-opt.fig')
    %}
end
if scheduling == 0
    brown_fcfs(d*24-(TS-1):d*24) = max(0,consumption_fcfs(d*24-(TS-1):d*24)+cooling_fcfs-r_p);
    brown_fcfs_s = sum(brown_fcfs);
    brown(d*24-(TS-1):d*24) = max(0,consumption(d*24-(TS-1):d*24)+cooling'-r_p);
    brown_s = sum(brown);
    compare1 = [sum(consumption_fcfs+cooling_fcfs)-sum(brown_fcfs), sum(consumption+cooling')-sum(brown);sum(brown_fcfs),sum(brown)]/1000*LN/360;
    figure;
    bar(compare1',0.3,'stacked')
    ylabel('energy consumption (KWh)');
    ylim([0,5]);
    hold on
    %plot(ones(1,2)*sum(r(TS*(day-1)+1:TS*day)),'-rs', 'LineWidth', 2)
    %title('Total power consumption and PV supply');
    legend('renewable','non-renewable');
    set(gca,'xticklabel',{'FCFS','Optimal'});
    %{
    compare2 = [(consumption)/sum(sum(b_record)),(consumption_flat)/sum(sum(b_flat_record))];
    compare2 = 1./compare2;
    figure;
    bar(compare2',0.3)
    ylim([0,1.2*max(compare2)]);
    title('Energy efficiency');
    ylabel('jobs/KWh(non-renewable)');
    %legend('grid power','renewable')
    set(gca,'xticklabel',{'Optimal','Simple'});
    
    figure;
    bar([sum(c_record(1:TS)+cooling(1:TS)'), sum(r(1:TS))],0.3)
    ylim([0,1200]);
    title('Consumption/supply');
    ylabel('KWh');
    %legend('total power consumption','total renewable supply')
    set(gca,'xticklabel',{'Consumption','Supply'});
    %}
    figure;
    bar([sum(r_p)-sum(consumption_fcfs)-sum(cooling)+brown_fcfs_s ,sum(r_p)-sum(consumption)-sum(cooling)+brown_s]/1000*LN/360,0.3)
    ylim([0,5]);
    %title('Surplus renewable');
    ylabel('KWh');
    %legend('total power consumption','total renewable supply')
    set(gca,'xticklabel',{'FCFS','Optimal'});
end
%{
    % baseline 2: EDF
    [Y, I] = sort(E);
    spare(d*24-(TS-1):d*24) = C*ones(1,24) - a(d*24-(TS-1):d*24);   
    for i = 1:1:N
        remain = B(I(i));
        current = S(I(i))+(d-1)*24-1;
        while remain > 0
            current = current +1;
            if current > E(I(i))+(d-1)*24
                !echo infeasible!
                I(i)
                remain
                remain = 0;
            end
            if remain <= spare(current)
                b_edf(current,I(i)) = remain;
                spare(current) = spare(current) - remain;
                remain = 0;
                break;
            else
                b_edf(current,I(i)) = spare(current);
                remain = remain - spare(current);
                spare(current) = 0;
            end
        end
    end    
    b_edf_record(d*24-(TS-1):d*24,:) = b_edf(d*24-(TS-1):d*24,:);
    c_edf_record(d*24-(TS-1):d*24) = sum(b_edf_record(d*24-(TS-1):d*24,:)') + a(d*24-(TS-1):d*24);
    edf(d) = p(d*24-(TS-1):d*24)*max(0,c_edf_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)');
    consumption_edf(d) = sum(max(0,c_edf_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)'));

    % baseline 3: FLAT
    for i = 1:1:N
        b_flat(S(i)+(d-1)*24:E(i)+(d-1)*24,i) = B(i)/(E(i)-S(i)+1);
    end    
    b_flat_record(d*24-(TS-1):d*24,:) = b_flat(d*24-(TS-1):d*24,:);
    c_flat_record(d*24-(TS-1):d*24) = sum(b_flat_record(d*24-(TS-1):d*24,:)') + a(d*24-(TS-1):d*24);
    flat(d) = p(d*24-(TS-1):d*24)*max(0,c_flat_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)');
    consumption_flat(d) = sum(max(0,c_flat_record(d*24-(TS-1):d*24)'-r(d*24-(TS-1):d*24)'));
    utilization_flat(d) = 1-sum(max(0,r(d*24-(TS-1):d*24)-c_flat_record(d*24-(TS-1):d*24)))/sum(r(d*24-(TS-1):d*24));
    
    % baseline 4: Renewable-supply-oblivious
    cvx_begin
        variables b(TS,N) c(TS)
        minimize( p(d*24-(TS-1):d*24)*c )
        %minimize( p(d*24-(TS-1):d*24)*max(0,c-r(d*24-(TS-1):d*24)') - R*sum(A(d*24-(TS-1):d*24,:).*b)' )
        subject to
            c - sum(b')' >= a(d*24-(TS-1):d*24)';
            sum(A(d*24-(TS-1):d*24,:).*b)' == B';
            %c-r(d*24-(TS-1):d*24)' <= C/2*ones(TS,1); % bandwidth constraint
            %c >= min(r(d*24-(TS-1):d*24)', C); % free renewable
            c <= C;
            b >= 0;
    cvx_end
    c_oblivious_record(d*24-(TS-1):d*24) = c;
    b_oblivious_record(d*24-(TS-1):d*24,:) = b;
    oblivious(d) = cvx_optval;
    consumption_oblivious(d) = sum(max(0,c-r(d*24-(TS-1):d*24)'));
    utilization_oblivious(d) = 1-sum(max(0,r(d*24-(TS-1):d*24)'-c))/sum(r(d*24-(TS-1):d*24));
end


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
day = 1;


figure;
plot(1:TS,C*ones(1,TS),'k:',1:TS,c_record(TS*(day-1)+1:TS*day),'k',1:TS,a(TS*(day-1)+1:TS*day),'b',1:TS,sum(b_record(TS*(day-1)+1:TS*day,:)'),'g',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'r','LineWidth', 2);
xlabel('hour');
ylabel('demand');
xlim([1,TS]);
ylim([0,2*C])
set(gca,'XTick',[0:6:TS]);
legend('total capacity','total demand','interactive','batch','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', 'results/demand-opt.eps');
saveas(gcf,'results/demand-opt.fig')



if scheduling == 1 || scheduling == 0
    
    figure;
    a_p = a./c_record.*consumption;
    b_p = consumption - a_p;
    total = [a_p;b_p];
    bar(total','stacked');
    hold on
    plot(1:TS, r_p(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    legend('interactive','batch','renewable supply')
    %{
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    hold on
    bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    xlabel('hour');
    ylabel('cpu shares');
    ylim([0,1.5*C]);
    legend('capacity','renewable supply','interactive','batch')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-opt.eps');
    saveas(gcf,'results/solution-opt.fig')
    %}
end

if scheduling == 2 || scheduling == 0
    total_fcfs = [a;b_fcfs_record'];
    figure;
    bar(total_fcfs(:,TS*(day-1)+1:TS*day)','stacked');
    hold on
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    xlabel('hour');
    ylabel('both workload');
    ylim([0,2*C]);
    legend('interactive','job 1', 'job 2','job 3','capacity','renewable')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-fcfs.eps');
    saveas(gcf,'results/solution-fcfs.fig')
end

if scheduling == 3 || scheduling == 0
    total_edf = [a;b_edf_record'];
    figure;
    bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    hold on
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    xlabel('hour');
    ylabel('both workload');
    ylim([0,2*C]);
    legend('interactive','job 1', 'job 2','job 3','capacity','renewable')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-edf.eps');
    saveas(gcf,'results/solution-edf.fig')
end

if scheduling == 4 || scheduling == 0
    total_flat = [a;b_flat_record'];
    figure;
    bar(total_flat(:,TS*(day-1)+1:TS*day)','stacked');
    hold on
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    xlabel('hour');
    ylabel('both workload');
    ylim([0,2*C]);
    legend('interactive','job 1', 'job 2','job 3','capacity','renewable')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-flat.eps');
    saveas(gcf,'results/solution-flat.fig')
end

if scheduling == 5 || scheduling == 0
    total_oblivious = [a;b_oblivious_record'];
    figure;
    bar(total_oblivious(:,TS*(day-1)+1:TS*day)','stacked');
    hold on
    plot(1:TS,C*ones(1,TS),'k:',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'-rs', 'LineWidth', 2)
    xlabel('hour');
    ylabel('both workload');
    ylim([0,2*C]);
    legend('interactive','job 1', 'job 2','job 3','capacity','renewable')
    xlim([1,TS]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', 'results/solution-oblivious.eps');
    saveas(gcf,'results/solution-oblivious.fig')
end

if scheduling == 0
figure;
bar([optimal(1)/optimal(1),fcfs(1)/optimal(1),flat(1)/optimal(1);consumption(1)/consumption(1),consumption_fcfs(1)/consumption(1),consumption_flat(1)/consumption(1);utilization(1)/utilization(1),utilization_fcfs(1)/utilization(1),utilization_flat(1)/utilization(1)]);ylim([0,3]);
legend('optimal','FCFS', 'FLAT');
set(gca,'yticklabel',{'0%', '50%', '100%', '150%', '200%', '250%', '300%'},'xticklabel',{'power cost', 'grid power usage', 'renewable utilization'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/comparisons.eps');
saveas(gcf,'results/comparisons.fig')

% timestamp | renewable supply | electricity price | 
% total grid power consumption (optimal) | total renewable consumption 
% (optimal) | total power consumption (optimal) | total power cost (optimal)
% | total grid power consumption (FCFS) | total renewable consumption
% (FCFS) | total power consumption (FCFS) | total power cost (FCFS)
% | total grid power consumption (FLAT) | total renewable consumption
% (FLAT) | total power consumption (FLAT) | total power cost (FLAT)
%{
csvwrite('results\problem2.csv',[(1:TS);r(1:TS);p_original(1:TS);max(0,c_record(1:TS)-r(1:TS));min(c_record(1:TS),r(1:TS));c_record(1:TS);p_original(1:TS).*max(0,c_record(1:TS)-r(1:TS))...
    ;max(0,c_fcfs_record(1:TS)-r(1:TS));min(c_fcfs_record(1:TS),r(1:TS));c_fcfs_record(1:TS);p_original(1:TS).*max(0,c_fcfs_record(1:TS)-r(1:TS))...
    ;max(0,c_flat_record(1:TS)-r(1:TS));min(c_flat_record(1:TS),r(1:TS));c_flat_record(1:TS);p_original(1:TS).*max(0,c_flat_record(1:TS)-r(1:TS))]',1,0);
%}

end

%}
%}