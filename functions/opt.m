function [cost_opt,capacity_opt] = opt(T,IP,PP,OP,CO2_grid,a,au,con,BS,A,S,E,PUE,solar,P,GP,RP,BP,plt,location,C_flat)
% number of servers used by interactive workloads
a_opt_server = ceil(a./(au'*ones(1,T))/PP);
a_opt_power = a_opt_server*(IP + (PP-IP)*au);
a_remaining = a_opt_server*PP - a;

cvx_begin
    variables PV C b1(size(BS,1),T) b2(size(BS,1),T) SB(T) G1(T)
    minimize (24*365/T*(RP*max(0,(sum(a_opt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB' - G1')' + GP*G1 - BP*SB)  + P(3)*PV + (P(1)+P(2))*C)
    subject to
        PV >= 0;
        PV <= 3*C;
        C*OP >= max(a_opt_server + sum(b2,1)/PP);
        sum(b1,1) <= a_remaining*con;
        b1 >= 0;
        b2 >= 0;
        sum(A.*b1,2) + sum(A.*b2,2) == BS;
        sum(b1,2) + sum(b2,2) == BS;
        SB >= 0;
        SB <= PV*solar';
        sum(SB) >= sum(G1);
        G1 >= 0;
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;
C_opt = C;

cost_opt(1) = C* PP * P(1);
cost_opt(2) = C* PP * mean(PUE-1) * P(2);
cost_opt(3) = PV * P(3) - 24*365/T*BP*SB;
cost_opt(4) = C* PP * P(4);
cost_opt(5) =24*365/T*(RP*max(0,(sum(a_opt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE - PV*solar + SB' - G1')' + GP*G1);

capacity_opt(1) = C_opt* PP;
capacity_opt(2) = C_opt* PP * mean(PUE-1);
capacity_opt(3) = PV;
capacity_opt(4) = C_opt* PP;

% CO2 emission
CO2_opt = CO2_grid*sum((sum(a_opt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*PUE- PV*solar)*24*365/T;


if plt >= 1
    figure;
    bar([sum(a_opt_power,1);sum(b1*(PP-IP)/PP,1) + sum(b2,1);(sum(a_opt_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*PP*ones(1,T),'k',1:T,PV*solar,'r', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C_flat*1.5]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
    
    if plt >= 2
        figure;
        bar([a_opt_server;sum(b2,1)/PP]','stacked')
        hold on;
        plot(1:T,C*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('server number');
        xlim([1,T]);
        ylim([0,C*1.1]);
    end
end