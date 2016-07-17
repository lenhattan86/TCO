function [cost_flat,capacity_flat] = no_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,P,GP,plt,location)
% number of servers used by interactive workloads
a_flat_server = ceil(a./(au'*ones(1,T))/PP);
a_flat_power = max(a_flat_server,[],2)*ones(1,T)*IP + a_flat_server*(PP-IP)*au;

% batch job scheduling: run with flat rate
b_flat = zeros(size(BS,1),T);
for i = 1:1:size(BS,1)
    b_flat(i,S(i):E(i)) = ones(1,E(i)-S(i)+1)*(BS(i)/(E(i)-S(i)+1));
end
% number of servers used by batch jobs
b_flat_server = sum(ceil(b_flat/PP),1);
b_flat_power = max(b_flat_server,[],2)*ones(1,T)*IP + b_flat_server*(PP-IP);

% decide the total number of server required by batch
C_flat = (sum(max(a_flat_server,[],2)) + sum(max(b_flat_server,[],2)))/OP;

cost_flat(1) = C_flat * PP * P(1);
cost_flat(2) = C_flat * PP * mean(PUE-1) * P(2);
cost_flat(3) = 0 * P(3);
cost_flat(4) = C_flat * PP * P(4);
cost_flat(5) = GP*((sum(a_flat_power,1) + sum(b_flat_power,1)).*PUE)'*24*365/T;

capacity_flat(1) = C_flat* PP;
capacity_flat(2) = C_flat* PP * mean(PUE-1);
capacity_flat(3) = 0;
capacity_flat(4) = C_flat* PP;

CO2_flat = CO2_grid*(sum(a_flat_power,1) + sum(b_flat_power,1))*PUE'*24*365/T;

if plt >= 1
    figure;
    bar([sum(a_flat_power,1);sum(b_flat_power,1);(sum(a_flat_power,1) + sum(b_flat_power,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C_flat*PP*ones(1,T),'k', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_flat*PP*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
    if plt >= 2
        figure;
        bar([a_flat_server;sum(b_flat_server,1)]','stacked')
        hold on;
        plot(1:T,C_flat*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C_flat*1.1]);
    end
end