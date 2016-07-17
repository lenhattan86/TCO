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
T = 24*7; % planning period length, in hours
IP = 0.1; % idle power consumption of a server, in kW
PP = 0.4; % peak power consumption of a server, in kW
CO2_grid = 0.5; % kg CO2/kWh of grid power consumption
DC_power = 500; % Data center power capacity, in kW
OP = 0.8; % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.8; % maximum utilization after consolidation

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
au = 0.3; % average utilization of interactive workloads
IPU = 0.9; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 3; % peak to mean ratio, should be smaller than current PMR
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
end

% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.2; % batch job ratio, compared with interactive workload
[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',11.99,23.99,'Uniform',1,2,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties

% uncontrolled renewable generation, e.g., PV, wind
RN = 1; % number of renewable generation
RL = {'traces\solar-one-week-2.csv'};
R = zeros(IN,T);
RP = 0.2; % capacity factor
RR = [0, 500]; % renewable capacity range, in kW
for i = 1:1:RN
    R(i,:) = trace_process(char(RL(i)), T, 12, 1, 1, RP(i));
end
RC = [3200/20, 0.01+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]

% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
CRC = [1000/20+SR*200/4, 0.02+0.12]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]

% cooling efficiency
PUE = 1.3*ones(1,T);
PUE_low = 1.3*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
GP = 0.12*ones(1,T); % grid electricity price
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          SOLVER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_option = 0;
[cost_houston,capacity_houston,G] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),RR,RC,SR,CRR,CRC,P,GP,RP,BP,plot_option,'results\houston');
[cost_houston_flat,capacity_houston_flat,G_flat] = houston_flat(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),RR,RC,SR,CRR,CRC,P,GP,RP,BP,plot_option,'results\houston_flat');
[cost_houston_noPV,capacity_houston_noPV] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),[0,0],RC,SR,CRR,CRC,P,GP,RP,BP,plot_option,'results\houston_noPV');
[cost_houston_givenPV,capacity_houston_givenPV] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),[250,250],RC,SR,CRR,CRC,P,GP,RP,BP,plot_option,'results\houston_givenPV');
[cost_houston_givenGE,capacity_houston_givenGE] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),RR,RC,SR,[1400,1400],CRC,P,GP,RP,BP,plot_option,'results\houston_givenGE');
[cost_houston_givenBoth,capacity_houston_givenBoth] = houston(T,IP,PP,OP,IT,CO2_grid,a,au,con,BS,A,S,E,PUE,R(1,:),[250,250],RC,SR,[1400,1400],CRC,P,GP,RP,BP,plot_option,'results\houston_givenBoth');

%{
C = IT;
solar = R(1,:);
location = 'results\houston';
plt = 2;
range = RR;


a_houston_server = ceil(a./(au'*ones(1,T))/PP);
a_houston_power = a_houston_server*(IP + (PP-IP)*au);
a_remaining = a_houston_server*PP*(con-au);
cvx_begin
    variables PV GE G(T) b1(size(BS,1),T) b2(size(BS,1),T)
    minimize (24*365/T*sum(RC(2)*PV*solar' + CRC(2)*G)  + RC(1)*PV + CRC(1)*GE)
    subject to
        PV >= range(1);
        PV <= range(2);
        C*OP/PP >= max(a_houston_server + sum(b2,1)/con/PP);
        sum(b1,1) <= sum(a_remaining,1);
        b1 >= 0;
        b2 >= 0;
        sum(A.*b1,2) + sum(A.*b2,2) == BS;
        sum(b1,2) + sum(b2,2) == BS;
        (sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(IP + (PP-IP)*con)).*PUE <= PV*solar + G';
        G >= 0;
        G <= GE*ones(T,1);
        GE >= 0;
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;
cvx_optval
PV
GE

cost_houston(1) = C * P(1); % IT install cost
cost_houston(2) = C * mean(PUE-1) * P(2); % cooling capaicty cost
cost_houston(3) = C * P(4); % server cost
cost_houston(4) = PV * RC(1); % PV install cost
cost_houston(5) = GE * CRC(1); % GE install cost
cost_houston(6) = 24*365/T*sum(RC(2)*PV*solar'); % PV O&M cost
cost_houston(7) = 24*365/T*sum(CRC(2)*G); % GE O&M cost

capacity_houston(1) = PV;
capacity_houston(2) = GE;

if plt >= 1
    figure;
    bar([sum(a_houston_power,1);sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(IP + (PP-IP)*con);(sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(IP + (PP-IP)*con)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*ones(1,T),'k',1:T,PV*solar+G','r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
    
    if plt >= 2
        figure;
        bar([a_houston_server;sum(b2,1)/PP/con]','stacked')
        hold on;
        plot(1:T,C/PP*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C/PP*1.1]);
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          OUTPUT                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
bar([cost_houston;cost_houston_flat;cost_houston_noPV;cost_houston_givenPV;cost_houston_givenGE;cost_houston_givenBoth],0.2,'stacked');
ylabel('annual expenditure ($)');
legend('PV installed cost','PV O&M cost','GE installed cost','GE O&M cost');
ylim([0,sum(cost_houston_givenBoth)*1.3]);
set(gca,'xticklabel',{'Our solution','no demand shaping','no PV','PV = 250kW','GE=1.4mW','PV = 250kW,GE=1.4mW'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\expenditure_comp_houston.eps');
saveas(gcf,'results\expenditure_comp_houston.fig');

figure;
bar([capacity_houston;capacity_houston_flat;capacity_houston_noPV;capacity_houston_givenPV;capacity_houston_givenGE;capacity_houston_givenBoth]);
ylabel('capacity (kW)');
legend('PV','GE');
ylim([0,max(capacity_houston_givenBoth)*1.3]);
set(gca,'xticklabel',{'Our solution','no demand shaping','no PV','PV = 250kW','GE=1.4mW','PV = 250kW,GE=1.4mW'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results\capacity_comp_houston.eps');
saveas(gcf,'results\capacity_comp_houston.fig');

save('results\houston.mat');


%{
figure;
bar([CO2_flat;CO2_PV],0.2);
ylabel('annual CO2 emission (kg)');
set(gca,'xticklabel',{'No integration','PV integration'});
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'CO2_comp.eps');
saveas(gcf,'CO2_comp.fig');
%}