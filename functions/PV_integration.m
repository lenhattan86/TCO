function [cost_PV,capacity_PV] = PV_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,solar,P,GP,BP,plt,location)
% number of servers used by interactive workloads
a_PV_server = ceil(a./(au'*ones(1,T))/PP);
a_PV_power = max(a_PV_server,[],2)*ones(1,T)*IP + a_PV_server*(PP-IP)*au;

% batch job scheduling: run with flat rate
b_PV = zeros(size(BS,1),T);
for i = 1:1:size(BS,1)
    b_PV(i,S(i):E(i)) = ones(1,E(i)-S(i)+1)*(BS(i)/(E(i)-S(i)+1));
end
% number of servers used by batch jobs
b_PV_server = sum(ceil(b_PV/PP),1);
b_PV_power = max(b_PV_server,[],2)*ones(1,T)*IP + b_PV_server*(PP-IP);

% decide the total number of server required by batch
C_PV = (sum(max(a_PV_server,[],2)) + sum(max(b_PV_server,[],2)))/OP;

% decide the PV capacity
PV = (sum(a_PV_power,1) + sum(b_PV_power,1))*PUE'/sum(solar);

cost_PV(1) = C_PV* PP * P(1);
cost_PV(2) = C_PV* PP * mean(PUE-1) * P(2);
cost_PV(3) = PV*P(3)- BP*(max(0,PV*solar-(sum(a_PV_power,1) + sum(b_PV_power,1)).*PUE ))' *24*365/T;
cost_PV(4) = C_PV* PP * P(4);
cost_PV(5) = GP*(max(0,(sum(a_PV_power,1) + sum(b_PV_power,1)).*PUE - PV*solar))' *24*365/T;

capacity_PV(1) = C_PV* PP;
capacity_PV(2) = C_PV* PP * mean(PUE-1);
capacity_PV(3) = PV;
capacity_PV(4) = C_PV* PP;

% CO2 emission
CO2_PV = CO2_grid*sum((sum(a_PV_power,1) + sum(b_PV_power,1)).*PUE- PV*solar)*24*365/T;

if plt >= 1
    figure;
    bar([sum(a_PV_power,1);sum(b_PV_power,1);(sum(a_PV_power,1) + sum(b_PV_power,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C_PV*PP*ones(1,T),'k',1:T,PV*solar,'r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_PV*PP*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
    
    if plt >= 2
        figure;
        bar([a_PV_server;sum(b_PV_server,1)]','stacked')
        hold on;
        plot(1:T,C_PV*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C_PV*1.1]);
    end
end