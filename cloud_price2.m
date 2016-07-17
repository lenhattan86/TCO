has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_clear;
clc;
clear;
%{
choice = menu('Pick one setting','Optimal','Night','FCFS','Flat','All');
switch choice
    case 1,
        scheduling = 1;
    case 2,
        scheduling = 2;
    case 3,
        scheduling = 3;
    case 4,
        scheduling = 4;
    case 5,
        scheduling = 0;
end
%}
scheduling = 0;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          INPUT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constant input
D = 1; % total days
TS = 24; % number of time slots in one day
T = TS*D; % total time slots of interests
TR = 24.6; % service demand of each request, in CPUshares*sec
PV = 518400; % peak PV hourly supply in CPUshares*sec
TN = 12; % number of slots per hour
SC = 30; % scale Peak PV to 120kW
% hourly temperature, from real traces
t_raw = load('traces\hourly-temperature.txt');
%COP = 6*ones(TS,1);

RAT = 35;
MaxBlowerCFM=10000;
MaxBlowerPower=11;
NB=1;
beta=3;
% hourly temperature
OAT = (t_raw - 32)/1.8;
y = 4*400;
COP= 4.246 + 4766*y./(OAT.^(1/2)*(y*y-1548*y+8.389e5)+5.863e5 );
AirCp = 1.006; % air capacity, kW/degC
MaxBlowerFlowRate = MaxBlowerCFM * 1.16/2119; % in kg/sec
% PUE
%{
PUEHP = load('traces\PUEHP.csv');
PUE_sum = zeros(1,24);
for d = 1:1:365
    for t = 1:1:24
        PUE_sum(t) = PUE_sum(t) + PUEHP((d-1)*24+t);
    end
end
PUE = PUE_sum/365;
%}
%PUE = [1.16,1.17,1.16,1.20,1.22,1.22,1.24,1.26,1.35,1.32,1.25,1.30,1.29,...
%    1.35,1.32,1.40,1.40,1.25,1.29,1.30,1.28,1.29,1.18,1.13];
%PUE = ones(1,24);

% PV solar traces
PAPV = load('traces\solar-one-week.csv'); % PV solar traces
for t = 1:1:T % average to each hour
    r(t) = max(0,sum(PAPV((t)*12-11:(t)*12)));
end
MR = 3600/SC;
r = r/max(r)*PV;

% transactional workload traces
%L = load('traces\VDR\VDR.txt'); % VDR workload trace for transactional workload
L = load('traces\SAPnew/sapTrace.tab');  % SAP workload trace for transactional workload

%{
figure;
plot(1/TN:1/TN:TS*7,L(1:TS*7*TN,4)/max(L(1:TS*7*TN,4)))
xlabel('hour');
ylabel('normalized CPU demand');
xlim([0,TS*7]);
ylim([0,1]);
set(gca,'XTick',[0:TS:TS*7]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 6.0 2.0]);
print ('-depsc', 'results\cpu_demand.eps');
%}


at = L(1:T*TN,4);

at = floor(at/max(at)*800)/10+0.1; % scale the average 
%{
fid = fopen('exp\run.sh','w');
fprintf(fid,'./httperf --server=192.168.199.159 --period=p');
interval = 3600/TN/SC;
for t = 1:1:T*TN
    fprintf(fid,'%.1f,%d',at(t),interval);
    if t < T*TN
        fprintf(fid,',');
    end
end
fprintf(fid,' --uri=/cgi-bin/sysbench.cgi --num-conn=%d --timeout=200 --print-reply',sum(at)*interval);
fclose(fid);
%}
for t = 1:1:T % average to each hour
    a(t) = mean(at(t*12-11:t*12))*24.6*3600/SC;
end


%{
figure;
plot(1:T,r,'b',1:T,a,'r');
xlabel('hour');
ylabel('power');
xlim([0,24]);
legend('renewable supply','transactional workload');
set(gca,'XTick',[0:6:24]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', 'results\power.eps');
saveas(gcf,'results\power.fig')
%}

% capacity cap, from real traces and model % grid electricity price, from real traces
% C = max(r(1:TS))/mean(PUE);
C = 100*PV/MR;
% grid electricity price, from real traces
p = load('traces\pa-electric-price.txt');
p = p/10; % change to cents/kWh
%p = mean(p)*ones(1,24);
%{
figure;
plot(1/24:1/24:7,p);
xlabel('day');
ylabel('electricity price (cents/kWh)');
xlim([0,7]);
ylim([0,7]);
set(gca,'XTick',[0:1:7]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 6.0 2.0]);
print ('-depsc', 'results\electricityprice.eps');
%saveas(gcf,'results\electricity-price.fig')
%}

p_original = p;
if min(p) <= 0
    p = p - min(p);
end

% batch workload
BN = 10; % average number of batch job arrivals per timeslot
BM = 0.5; % batch job ratio, compared with interactive workload
N = BN*T; % total number of batch jobs
[A,B,S,E] = batch_job_generator(T,N,'Uniform',11.99,11.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)));
A = A';
B = B';
NB = PV/MR;

% price for both applications
size_per_instance = 2;
p_opt = 8*ones(TS,1);
p_fcfs = 12*ones(TS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          BEST EFFORT                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spare = zeros(1,TS*D);
b_fcfs = zeros(TS*D,N);
b_fcfs_record = zeros(TS*D,N);
c_fcfs_record = zeros(1,TS*D);

for d = 1:1:D
    if scheduling == 3 || scheduling == 0
        % baseline 3: FCFS
        [Y, I] = sort(S);
        spare(d*24-23:d*24) = C*ones(1,24) - a(d*24-23:d*24);   
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
                    b_fcfs(current,I(i)) = remain;
                    spare(current) = spare(current) - remain;
                    remain = 0;
                    break;
                else
                    b_fcfs(current,I(i)) = spare(current);
                    remain = remain - spare(current);
                    spare(current) = 0;
                end
            end
        end    
        b_fcfs_record(d*24-23:d*24,:) = b_fcfs(d*24-23:d*24,:);
        c_fcfs_record(d*24-23:d*24) = sum(b_fcfs_record(d*24-23:d*24,:)') + a(d*24-23:d*24);
        cvx_begin
            variables s(TS) % b is per CPU share
            minimize( p(d*TS-(TS-1):d*TS)*max(0,c_fcfs_record(d*TS-(TS-1):d*TS)'+(max(0,c_fcfs_record(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta)) - r(d*TS-(TS-1):d*TS)'))
            subject to
                s >= 0;
                s <= 1;
        cvx_end
        cooling_fcfs(d*24-23:d*24) = max(0,c_fcfs_record(d*TS-(TS-1):d*TS)' - NB*AirCp*MaxBlowerFlowRate * s.*max(0,RAT' - OAT'))./(COP') + NB*MaxBlowerPower * pow_pos(s,beta);
        fcfs(d) = p(d*24-23:d*24)*max(0,c_fcfs_record(d*24-23:d*24)'+ cooling_fcfs(d*24-23:d*24)'-r(d*24-23:d*24)');
        consumption_fcfs(d) = sum(max(0,c_fcfs_record(d*24-23:d*24)'+ cooling_fcfs(d*24-23:d*24)'-r(d*24-23:d*24)'));
        cost_fcfs(d) = sum(p(d*24-23:d*24)'.*max(0,c_fcfs_record(d*24-23:d*24)'+ cooling_fcfs(d*24-23:d*24)'-r(d*24-23:d*24)'));
        cost_other(d) = 2*cost_fcfs(d);
        revenue_fcfs(d) = sum(A(d*24-23:d*24,:).*b_fcfs,2)'/size_per_instance*p_fcfs;
        profit_fcfs(d) = revenue_fcfs(d)-cost_fcfs(d)-cost_other(d);
        utilization_fcfs(d) = sum(r(d*24-23:d*24))-sum(max(0,r(d*24-23:d*24)-c_fcfs_record(d*24-23:d*24)- cooling_fcfs(d*24-23:d*24)));
    end
end
total_cost_fcfs = sum(cost_fcfs);
total_cost_other = sum(cost_other);
total_revenue_fcfs = sum(revenue_fcfs);
total_profit_fcfs = sum(profit_fcfs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_record = zeros(T,N);
for d = 1:1:D
    cvx_begin
        variables b(TS,N) c(TS) s(TS)
        minimize( p(d*24-23:d*24)*max(0,a(d*24-23:d*24)+sum(b')+...
            max(0,a(d*TS-(TS-1):d*TS) + sum(b') - NB*AirCp*MaxBlowerFlowRate * s'.*max(0,RAT - OAT))./(COP) + NB*MaxBlowerPower * pow_pos(s',beta)...
                -r(d*24-23:d*24)  )' ) 
        %minimize( p(d*24-23:d*24)*max(0,c-r(d*24-23:d*24)') - R*sum(A(d*24-23:d*24,:).*b)' )
        subject to
            %c - sum(b')'.*PUE' >= a(d*24-23:d*24)'.*PUE';
            sum(A(d*24-23:d*24,:).*b)' == mean(p_fcfs)/mean(p_opt)*B';
            sum(b)' == mean(p_fcfs)/mean(p_opt)*B';
            %c-r(d*24-23:d*24)' <= C/2*ones(TS,1); % bandwidth constraint
            %c >= min(r(d*24-23:d*24)', C); % free renewable
            a(d*24-23:d*24)+sum(b') <= C;
            b >= 0;
            s >= 0;
            s <= 1;
    cvx_end
    if strcmp(cvx_status, 'Solved') == 0
        cvx_status
    end
    s_record = s;
    b_record(d*24-23:d*24,:) = b;
    cooling(d*24-23:d*24) = max(0,a(d*TS-(TS-1):d*TS) + sum(b') - NB*AirCp*MaxBlowerFlowRate * s'.*max(0,RAT - OAT))./(COP) + NB*MaxBlowerPower * pow_pos(s',beta);
    c_record(d*24-23:d*24) = a(d*24-23:d*24)+sum(b')+cooling(d*24-23:d*24);
    optimal(d) = cvx_optval;
    cost(d) = cvx_optval;
    revenue(d) = sum(A(d*24-23:d*24,:).*b,2)'/size_per_instance*p_opt;
    profit(d) = revenue(d)-cost(d)-cost_other(d);
    consumption(d) = sum(max(0,c_record(d*24-23:d*24)'-r(d*24-23:d*24)'));
    utilization(d) = sum(r(d*24-23:d*24))-sum(max(0,r(d*24-23:d*24)'-c_record(d*24-23:d*24)'));
end
total_cost = sum(cost);
total_revenue = sum(revenue);
total_profit = sum(profit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

day = 1;

%{
figure;
plot(1:TS,C*ones(1,TS),'k:',1:TS,c_record(TS*(day-1)+1:TS*day),'k',1:TS,a(TS*(day-1)+1:TS*day),'b',1:TS,sum(b_record(TS*(day-1)+1:TS*day,:)'),'g',1:TS,r_unchanged(TS*(day-1)+1:TS*day),'r','LineWidth', 2);
xlabel('hour');
ylabel('demand');
xlim([1,TS]);
ylim([0,2*C])
set(gca,'XTick',[0:6:TS]);
legend('total capacity','total demand','interactive','batch','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', 'results\demand-opt.eps');
saveas(gcf,'results\demand-opt.fig')
%}


if scheduling == 1 || scheduling == 0
    total = [a;sum(b_record');cooling]/PV*MR;
    figure;    
    h1 = bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    %title('Optimal')
    ylim([0,2*C/PV*MR]);
    xlabel('hour');
    ylabel('power (kW)');
    xlim([0,TS+1]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    hold on
    ah1 = gca;
    l1 = legend(ah1,h1,'reserved instances','flexible instances','cooling',2);
    h2 = plot(0:TS+1,C*ones(1,TS+2)/PV*MR,'k:',1:TS,r(TS*(day-1)+1:TS*day)/PV*MR,'-r', 'LineWidth', 2);
    ah2 = gca;
    ah2=axes('position',get(gca,'position'), 'visible','off');
    l2 = legend(ah2,h2,'IT capacity','renewable',1);
    LEG = findobj(l1,'type','text');
    set(LEG,'FontSize',8)
    LEG = findobj(l2,'type','text');
    set(LEG,'FontSize',8)
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\cloud2-opt.eps');
    saveas(gcf,'results\cloud2-opt.fig')
end




if scheduling == 3 || scheduling == 0
    total = [a;sum(b_fcfs_record');cooling_fcfs]/PV*MR;
    figure;
    h1 = bar(total(:,TS*(day-1)+1:TS*day)','stacked');
    %title('Optimal')
    ylim([0,2*C/PV*MR]);
    xlabel('hour');
    ylabel('power (kW)');
    xlim([0,TS+1]);
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    hold on
    ah1 = gca;
    l1 = legend(ah1,h1,'reserved instances','flexible instances','cooling',2);
    h2 = plot(0:TS+1,C*ones(1,TS+2)/PV*MR,'k:',1:TS,r(TS*(day-1)+1:TS*day)/PV*MR,'-r', 'LineWidth', 2);
    ah2 = gca;
    ah2=axes('position',get(gca,'position'), 'visible','off');
    l2 = legend(ah2,h2,'IT capacity','renewable',1);
    LEG = findobj(l1,'type','text');
    set(LEG,'FontSize',8)
    LEG = findobj(l2,'type','text');
    set(LEG,'FontSize',8)
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\cloud2-be.eps');
    saveas(gcf,'results\cloud2-be.fig')
end

figure;
bar([total_cost,total_cost_other,total_revenue,total_profit;total_cost_fcfs,total_cost_other,total_revenue_fcfs,total_profit_fcfs]);
legend('electricity cost','other cost','revenue','profit');
set(gca,'xticklabel',{'Our scheme','Best Effort'});
saveas(gcf,'results\cloud2-comp.fig')
%{
if scheduling == 0
    [consumption(1),consumption_night(1),consumption_fcfs(1),consumption_flat(1)]/PV*MR
    1-consumption(1)./[consumption(1),consumption_night(1),consumption_fcfs(1),consumption_flat(1)]
    [consumption(1)+utilization(1),consumption_night(1)+utilization_night(1),consumption_fcfs(1)+utilization_fcfs(1),consumption_flat(1)+utilization_flat(1)]/PV*MR
    (consumption(1)+utilization(1))./[consumption(1)+utilization(1),consumption_night(1)+utilization_night(1),consumption_fcfs(1)+utilization_fcfs(1),consumption_flat(1)+utilization_flat(1)]-1
    %{
    figure;
    x1 = [0.9,1.9,2.9,3.9];
    h1 = bar(x1,[consumption(1),consumption_night(1),consumption_fcfs(1),consumption_flat(1);utilization(1),utilization_night(1),utilization_fcfs(1),utilization_flat(1)]'/PV*MR,0.2,'stacked');
    %set(h1(1),'facecolor','g'); 
    %set(h1(2),'facecolor','r');
    hold on
    x2 = [1.1,2.1,3.1,4.1];
    h2 = bar(x2,[sum(a(1:TS)),sum(a(1:TS)),sum(a(1:TS)),sum(a(1:TS)); sum(sum(b_record(1:TS,:))),sum(sum(b_night_record(1:TS,:))),sum(sum(b_fcfs_record(1:TS,:))),sum(sum(b_flat(1:TS,:)));sum(cooling(1:TS)),sum(cooling_night(1:TS)),sum(cooling_fcfs(1:TS)),sum(cooling_flat(1:TS))]'/PV*MR,0.2,'stacked');
    set(h2(1),'facecolor','b'); 
    set(h2(2),'facecolor','c');
    set(h2(3),'facecolor','y');
    %title('Power consumption comparison')
    ylim([0,2000]);
    ylabel('power (kW)');
    set(gca,'XTick',[1:1:4], 'FontSize', 8);
    set(gca,'xticklabel',{'Optimal','Night', 'BE','Flat'});
    ah1 = gca;
    l1 = legend(ah1,h1,'grid power','renewable',2);
    ah2=axes('position',get(gca,'position'), 'visible','off');
    l2 = legend(ah2,h2,'interactive','batch job', 'cooling power',1);
    LEG = findobj(l1,'type','text');
    set(LEG,'FontSize',8)
    LEG = findobj(l2,'type','text');
    set(LEG,'FontSize',8)
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\a1-comp1.eps');
    print(gcf, '-dpdf', 'results\a1-comp2.pdf');
    %saveas(gcf,'results\a1-comp1.fig')
    %}
    figure;
    x1 = [0.9,1.9,2.9,3.9];
    h1 = bar([consumption(1),consumption_night(1),consumption_fcfs(1),consumption_flat(1);utilization(1),utilization_night(1),utilization_fcfs(1),utilization_flat(1)]'/PV*MR,0.4,'stacked');
    set(h1(1),'facecolor','b'); 
    set(h1(2),'facecolor','g');
    hold on
    ylim([0,2000]);
    ylabel('power (kWh)');
    set(gca,'XTick',[1:1:4], 'FontSize', 8);
    set(gca,'xticklabel',{'Optimal','Night', 'BE','Flat'});
    legend('grid power','renewable');
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\a1-comp1.eps');
    
    %saveas(gcf,'results\a1-comp1.fig')
    %figure;
    %bar(5000./([sum(B)/consumption(1),sum(B)/consumption_night(1),sum(B)/consumption_fcfs(1),sum(B)/consumption_flat(1)]'*PV/MR),0.3,'stacked');
    5000./([sum(B)/consumption(1),sum(B)/consumption_night(1),sum(B)/consumption_fcfs(1),sum(B)/consumption_flat(1)]'*PV/MR)
%{    
    %title('Efficiency comparison')
    %bar([utilization(1),consumption(1),sum(cooling(1:24));utilization_ncooling(1),consumption_ncooling(1),sum(cooling_ncooling(1:24));utilization_fcfs(1),consumption_fcfs(1),sum(cooling_fcfs(1:24));utilization_flat(1),consumption_flat(1),sum(cooling_flat(1:24))],'stacked');
    %ylim([0,8e6/PV*MR]);
    %legend('renewable','grid power');
    ylabel('kWh grid power/job unit');
    %set(gca,'yticklabel',{'0%', '50%', '100%', '150%', '200%', '250%', '300%'},'xticklabel',{'power cost', 'power usage', 'renewable utilization'});
    set(gca,'xticklabel',{'Optimal','Night', 'BE','Flat'});
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\a1-comp2.eps');
    print(gcf, '-dpdf', 'results\a1-comp2.pdf');
    %saveas(gcf,'results\a1-comp2.fig')

    
    figure;
    bar([cost(1),cost_night(1),cost_fcfs(1),cost_flat(1)]/100/PV*MR,0.3);
    1-cost(1)./[cost(1),cost_night(1),cost_fcfs(1),cost_flat(1)]
    %title('Efficiency comparison')
    %bar([utilization(1),consumption(1),sum(cooling(1:24));utilization_ncooling(1),consumption_ncooling(1),sum(cooling_ncooling(1:24));utilization_fcfs(1),consumption_fcfs(1),sum(cooling_fcfs(1:24));utilization_flat(1),consumption_flat(1),sum(cooling_flat(1:24))],'stacked');
    %ylim([0,8e6/PV*MR]);
    %legend('renewable','grid power');
    ylabel('energy cost ($)');
    %set(gca,'yticklabel',{'0%', '50%', '100%', '150%', '200%', '250%', '300%'},'xticklabel',{'power cost', 'power usage', 'renewable utilization'});
    set(gca,'xticklabel',{'Optimal','Night', 'FCFS','FLAT'});
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.2]);
    print ('-depsc', 'results\a1-comp3.eps');
%}    
    figure;
    comp3(1,:) = [cost(1),cost_night(1),cost_fcfs(1),cost_flat(1)]./max([cost(1),cost_night(1),cost_fcfs(1),cost_flat(1)]);
    comp3(2,:) = (1./([sum(B)/consumption(1),sum(B)/consumption_night(1),sum(B)/consumption_fcfs(1),sum(B)/consumption_flat(1)]))./max(1./([sum(B)/consumption(1),sum(B)/consumption_night(1),sum(B)/consumption_fcfs(1),sum(B)/consumption_flat(1)]));
    bar(comp3);
    1-cost(1)./[cost(1),cost_night(1),cost_fcfs(1),cost_flat(1)]
    ylim([0,2]);
    %title('Efficiency comparison')
    %bar([utilization(1),consumption(1),sum(cooling(1:24));utilization_ncooling(1),consumption_ncooling(1),sum(cooling_ncooling(1:24));utilization_fcfs(1),consumption_fcfs(1),sum(cooling_fcfs(1:24));utilization_flat(1),consumption_flat(1),sum(cooling_flat(1:24))],'stacked');
    %ylim([0,8e6/PV*MR]);
    legend('Optimal','Night', 'BE','FLAT');
    ylabel('normalized value');
    %set(gca,'yticklabel',{'0%', '50%', '100%', '150%', '200%', '250%', '300%'},'xticklabel',{'power cost', 'power usage', 'renewable utilization'});
    set(gca,'xticklabel',{'energy cost','CO2 emission'});
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.1]);
    print ('-depsc', 'results\a1-comp4.eps');

end
%}










