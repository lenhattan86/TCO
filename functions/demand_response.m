function [b,new_peak,power] = demand_response(nz,colocate,ds,off,T,IP,PP,OP,C,a,au,con,BS,A,bu,PUE,solar,range1,RC,RE,range2,CRC,CRE,GP,p_peak,peak,p_cp,t_cp,has_warning,WP,plt,location)
    idle_power = C/PP*IP*ones(1,T);
    a_houston_server = ceil(a./(au'*ones(1,T))/PP);
    a_houston_power = a_houston_server*((PP-IP)*au);
    if colocate == 1
        a_remaining = a_houston_server*PP*(con-au);
    else
        a_remaining = zeros(1,T);
    end
    b_flat = BS./sum(A,2)*ones(1,T).*A;
    if has_warning == 0
        GP = GP + WP*p_cp*T/24/30;
    end
    cvx_begin
        variables PV GE G(T) b1(size(BS,1),T) b2(size(BS,1),T) D(T)
        minimize (24*30/T*sum(RC(2)*PV*solar' + CRC(2)*G + GP'.*D)  + p_peak*max(0,max(D)-peak) + p_cp*D(t_cp)*has_warning)
        subject to
            PV >= range1(1);
            PV <= range1(2);
            C/PP >= max(a_houston_server + sum(b2,1)/bu/PP);
            sum(b1,1) <= sum(a_remaining,1);
            b1 >= 0;
            b2 >= 0;
            if ds == 0
                b1 + b2 == b_flat;
            end
            if off == 1
                idle_power = (a_houston_server + sum(b2,1)/bu/PP)*IP;
            end
            sum(A.*b1,2) + sum(A.*b2,2) == BS;
            sum(b1,2) + sum(b2,2) == BS;
            (idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE <= D' + PV*solar + G';
            D >= 0;
            D <= C/OP;
            (idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE <= C/OP;
            G >= 0;
            G <= GE*ones(T,1);
            GE >= range2(1);
            GE <= range2(2);
            if nz == 1
                sum((idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE) <= sum(PV*solar + G');
            end
    cvx_end

    if strcmp(cvx_status, 'Solved') == 0
        cvx_status
    end
    b = b1+b2;
    new_peak = max(peak,max(D));
    power = zeros(7,T);
    power = [idle_power;sum(a_houston_power,1);sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP;(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*(PUE-ones(1,T));PV*solar;G';D'];
    
    if plt >= 1
        figure;
        bar([idle_power;sum(a_houston_power,1);sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1);(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*(PUE-ones(1,T))]','stacked')
        hold on;
        plot(1:T,C*ones(1,T),'k',1:T,PV*solar+G','r', 'LineWidth', 2)
        xlabel('hour');
        ylabel('kW');
        xlim([1,T]);
        ylim([0,C*1.5]);
        legend('idle power','delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
        %set(gca,'XTick',[0:6:TS], 'FontSize', 8);
        set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
        print ('-depsc', strcat(location,'.eps'));
        saveas(gcf,strcat(location,'.fig'))
    end
end

