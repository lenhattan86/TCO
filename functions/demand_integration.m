function [cost_demand,capacity_demand] = demand_integration(T,IP,PP,OP,CO2_grid,a,au,con,BS,A,S,E,PUE,solar,P,GP,RP,BP,plt,location,capacity)
% number of servers used by interactive workloads
a_demand_server = ceil(a./(au'*ones(1,T))/PP);
a_demand_power = a_demand_server*IP + a_demand_server*(PP-IP)*au;
a_remaining = a_demand_server*PP - a;

% batch job scheduling: run with flat rate
b_demand = zeros(size(BS,1),T);
for i = 1:1:size(BS,1)
    b_demand(i,S(i):E(i)) = ones(1,E(i)-S(i)+1)*(BS(i)/(E(i)-S(i)+1));
end
% number of servers used by batch jobs
b_demand_server = sum(ceil(b_demand/PP),1);
b_demand_power = max(b_demand_server,[],2)*ones(1,T)*IP + b_demand_server*(PP-IP);

% decide the total number of server required by batch
C_demand = capacity(4);
PV = capacity(3);
C = C_demand;

cvx_begin
    variables b1(size(BS,1),T) b2(size(BS,1),T) SB(T)
    minimize (24*365/T*(GP*max(0,(sum(a_demand_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB')' - BP*SB))
    subject to
        sum(b1,1) <= a_remaining*con;
        b1 >= 0;
        b2 >= 0;
        sum(A.*b1,2) + sum(A.*b2,2) == BS;
        sum(b1,2) + sum(b2,2) == BS;
        SB >= 0;
        SB <= PV*solar';
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;

cost_demand(1) = capacity(1) * P(1);
cost_demand(2) = capacity(2) * P(2);
cost_demand(3) = PV * P(3) - 24*365/T*BP*SB;
cost_demand(4) = C * P(4);
cost_demand(5) = 24*365/T*GP*max(0,(sum(a_demand_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB')';

capacity_demand = capacity;

% CO2 emission
CO2_demand = CO2_grid*sum((sum(a_demand_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE- PV*solar)*24*365/T;

if plt >= 1
    figure;
    bar([sum(a_demand_power,1);sum(b1*(PP-IP)/PP,1) + sum(b2,1);(sum(a_demand_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*ones(1,T),'k',1:T,PV*solar,'r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_demand*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))

    if plt >= 2
        figure;
        bar([a_demand_server;sum(b2,1)/PP]','stacked')
        hold on;
        plot(1:T,C*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C_demand*1.1]);
    end

end