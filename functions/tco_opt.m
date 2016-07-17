function [cost,b,EU,status,C,EC,RC] = tco_opt(T,a,BS,A,R,PUE,P_IT,P_ES,P_R,Loss,plt,location)

%%% TODO: consider charge and discharge independently

cvx_begin
    variables C EC ES0 RC(size(R,1)) b(size(BS,1),T) ES(T) EU(T)
    % C: IT capacity, in kW
    % EC: energy storage capacity, in kWh
    % ES0: initial energy storage level, in kWh
    % RC: renewable capacity, in kW
    % b: batch jobs
    % ES: energy storage level
    % Ch: charge of energy storage
    % Di: discharge of energy storage
    minimize( P_IT*C + P_ES*EC + P_R*RC )
    subject to
        (sum(a,1) + sum(b,1)).*PUE <= RC'*R - EU';
        sum(a,1) + sum(b,1) <= C*ones(1,T);
        C >= 0;
        EC >= 0;
        ES <= EU + Loss(1)*[ES0;ES(1:T-1)];
        ES >= 0;
        ES <= EC;
        ES(T) >= ES0
        b >= 0;
        sum(A.*b, 2) == BS;
        sum(b, 2) == BS;
cvx_end

if strcmp(cvx_status, 'Infeasible') == 1
    'Infeasible'
end

status = cvx_status
cost = cvx_optval;

if plt == 1
    figure;
    bar([sum(a,1);sum(b,1);(sum(a,1) + sum(b,1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*ones(1,T),'k',1:T,RC'*R,'r', 1:T,RC'*R-EU','y', 'LineWidth', 2)
    xlabel('hour');
    ylabel('kW');
    xlim([0,T]);
    ylim([0,max([C;RC])*1.2]);
    legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable','renewable+storage')
    %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
    print ('-depsc', strcat(location,'.eps'));
    saveas(gcf,strcat(location,'.fig'))
end


