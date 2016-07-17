function [cost_supply,capacity_supply] = supply_integration(T,IP,PP,OP,CO2_grid,a,au,BS,A,S,E,PUE,solar,P,GP,RP,BP,plt,location)
% number of servers used by interactive workloads
a_supply_server = ceil(a./(au'*ones(1,T))/PP);
a_supply_power = max(a_supply_server,[],2)*ones(1,T)*IP + a_supply_server*(PP-IP)*au;

% batch job scheduling: run with flat rate
b_supply = zeros(size(BS,1),T);
for i = 1:1:size(BS,1)
    b_supply(i,S(i):E(i)) = ones(1,E(i)-S(i)+1)*(BS(i)/(E(i)-S(i)+1));
end
% number of servers used by batch jobs
b_supply_server = sum(ceil(b_supply/PP),1);
b_supply_power = max(b_supply_server,[],2)*ones(1,T)*IP + b_supply_server*(PP-IP);

% decide the total number of server required by batch
C_supply = (sum(max(a_supply_server,[],2)) + sum(max(b_supply_server,[],2)))/OP;

% optimization for PV capacity and grid power usage
        %- BP*(max(0,PV*solar-(sum(a_supply_power,1) + sum(b_supply_power,1)).*PUE ))'...
cvx_begin
    variables PV SB(T) G1(T)
    minimize (24*365/T*(RP*max(0,(sum(a_supply_power,1) + sum(b_supply_power,1)).*PUE - PV*solar + SB'-G1')' + GP*G1 - BP*SB)  + P(3)*PV)
    subject to
        PV >= 0;
        PV <= 3*C_supply;
        SB >= 0;
        SB <= PV*solar';
        sum(SB) >= sum(G1);
        G1 >= 0;
cvx_end
if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;

cost_supply(1) = C_supply* PP * P(1);
cost_supply(2) = C_supply* PP * mean(PUE-1) * P(2);
cost_supply(3) = PV * P(3) - BP*SB*24*365/T;
cost_supply(4) = C_supply* PP * P(4);
cost_supply(5) = (RP*((sum(a_supply_power,1) + sum(b_supply_power,1)).*PUE - PV*solar + SB'-G1')'+ GP*G1)*24*365/T;

capacity_supply(1) = C_supply* PP;
capacity_supply(2) = C_supply* PP * mean(PUE-1);
capacity_supply(3) = PV;
capacity_supply(4) = C_supply* PP;


% CO2 emission
CO2_supply = CO2_grid*sum((sum(a_supply_power,1) + sum(b_supply_power,1)).*PUE- PV*solar)*24*365/T;

if plt == 1
    figure;
    bar([sum(a_supply_power,1);sum(b_supply_power,1);(sum(a_supply_power,1) + sum(b_supply_power,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C_supply*PP*ones(1,T),'k',1:T,PV*solar,'r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_supply*PP*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
end