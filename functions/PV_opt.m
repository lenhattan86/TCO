function [cost_PVopt,CO2_PVopt,C_PVopt,status] = PV_opt(T,IP,PP,CO2_grid,a,au,con,BS,A,S,E,PUE,solar,P,GP,BP,plt,location,C_flat)
% number of servers used by interactive workloads
a_PVopt_server = ceil(a./(au'*ones(1,T))/PP);
a_PVopt_power = a_PVopt_server*(IP + (PP-IP)*au);
a_remaining = a_PVopt_server*PP - a;



cvx_begin
    variables PV C b1(size(BS,1),T) b2(size(BS,1),T) SB(T)
    minimize (24*365/T*(GP*max(0,(sum(a_PVopt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB')' - BP*SB)  + P(3)*PV + (P(1)+P(2))*C)
    subject to
        PV >= 0;
        PV <= 3*C;
        C >= max(a_PVopt_server + sum(b2,1)/PP);
        sum(b1,1) <= a_remaining*con;
        b1 >= 0;
        b2 >= 0;
        sum(A.*b1,2) + sum(A.*b2,2) == BS;
        sum(b1,2) + sum(b2,2) == BS;
        SB >= 0;
        SB <= PV*solar';
        sum((sum(a_PVopt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE) <= sum(PV*solar);
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;
C_PVopt = C;

cost_PVopt(1) = C * P(1);
cost_PVopt(2) = C * P(2);
cost_PVopt(3) = 24*365/T*GP*max(0,(sum(a_PVopt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB')';
cost_PVopt(4) = PV * P(3) - 24*365/T*BP*SB;

% CO2 emission
CO2_PVopt = CO2_grid*sum((sum(a_PVopt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE- PV*solar)*24*365/T;


if plt >= 1
    figure;
    bar([sum(a_PVopt_power,1);sum(b1*(PP-IP)/PP,1) + sum(b2,1);(sum(a_PVopt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*PP*ones(1,T),'k',1:T,PV*solar,'r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_flat*PP*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))

    if plt >= 2
        figure;
        bar([a_PVopt_server;sum(b2,1)/PP]','stacked')
        hold on;
        plot(1:T,C*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C_flat*1.1]);
    end

end