function [cost,capacity,G,status] = houston(T,IP,PP,OP,C,CO2_grid,a,au,con,BS,A,S,E,PUE,solar,range1,RC,SR,range2,CRC,P,GP,RP,BP,plt,location)
a_houston_server = ceil(a./(au'*ones(1,T))/PP);
a_houston_power = a_houston_server*(IP + (PP-IP)*au);
a_remaining = a_houston_server*PP*(con-au);
cvx_begin
    variables PV GE G(T) b1(size(BS,1),T) b2(size(BS,1),T)
    minimize (24*365/T*sum(RC(2)*PV*solar' + CRC(2)*G)  + RC(1)*PV + CRC(1)*GE)
    subject to
        PV >= range1(1);
        PV <= range1(2);
        C*OP/PP >= max(a_houston_server + sum(b2,1)/con/PP);
        sum(b1,1) <= sum(a_remaining,1);
        b1 >= 0;
        b2 >= 0;
        sum(A.*b1,2) + sum(A.*b2,2) == BS;
        sum(b1,2) + sum(b2,2) == BS;
        (sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(IP + (PP-IP)*con)/con/PP).*PUE <= PV*solar + G';
        G >= 0;
        G <= GE*ones(T,1);
        GE >= range2(1);
        GE <= range2(2);
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
end
status = cvx_status;
cvx_optval

cost = [PV * RC(1), 24*365/T*sum(RC(2)*PV*solar'), GE * CRC(1), 24*365/T*sum(CRC(2)*G)];
capacity = [PV, GE];
%{
cost_houston(1) = C * P(1); % IT install cost
cost_houston(2) = C * mean(PUE-1) * P(2); % cooling capaicty cost
cost_houston(3) = C * P(4); % server cost
cost_houston(4) = PV * RC(1); % PV install cost
cost_houston(5) = GE * CRC(1); % GE install cost
cost_houston(6) = 24*365/T*sum(RC(2)*PV*solar'); % PV O&M cost
cost_houston(7) = 24*365/T*sum(CRC(2)*G); % GE O&M cost

capacity_houston(1) = PV;
capacity_houston(2) = GE;
%}

if plt >= 1
    figure;
    bar([sum(a_houston_power,1);sum(b1/PP*(PP-IP),1)+ sum(b2,1)/con/PP*(IP + (PP-IP)*con);(sum(a_houston_power,1) + sum(b1/PP*(PP-IP),1)+ sum(b2,1)/con/PP*(IP + (PP-IP)*con)).*(PUE-ones(1,T))]','stacked')
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
    
    
    if plt >= 3
        figure;
        bar([G']','stacked')
        hold on;
        plot(1:T,GE*ones(1,T),'k', 'LineWidth', 2)
        xlabel('hour');
        ylabel('kW');
        legend('Gas Engine generation','Gas Engine capacity')
        xlim([1,T]);
        ylim([0,C*1.5]);
    end
end

