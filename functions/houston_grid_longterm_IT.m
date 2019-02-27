function [cost,capacity,emission, total_demand, total_supply] = houston_grid_longterm_IT(colocate,nz,ds,off, ...
    T, N_y, ...
    IP,PP,OP, ...
    ITR,ITC, ...
    CO2_grid, a, au, ...
    con,BS,A,S,E,bu,PUE,solar, ...
    range1,RC,RE,SR,range2,CRC,CRE,P,GP,RP,BP,plt,location)

% idle_power = C/PP*IP*ones(T,N_y);
a_houston_server = ceil(a./(au'*ones(T,N_y))/PP); % a_houston_server T,N_y
a_houston_power = a_houston_server*((PP-IP)*au);
if colocate == 1
    a_remaining = a_houston_server*PP*(con-au);
else
    a_remaining = zeros(T,N_y);
end
for y=1:N_y
    b_flat(:,:,y) = BS(:,y)./sum(A(:,:,y),2)*ones(1, T).*A(:,:,y); % todo: repeat A -> (BN, T, N_y); % 
end

if and(range2(1) > 0 ,range2(1)==range2(2) )
    D_cap = 0;
end
% PV: PV capacity
% GE: Energy storage capacity
% D: Demand capacity
% C: IT capacity
cvx_begin
    variables PV(N_y) GE(N_y) G(T,N_y) b1(size(BS,1),T,N_y) b2(size(BS,1),T,N_y) D(T,N_y) C(N_y)
    minimize (sum(24*365/T*sum((ones(T,1)*RC(2,:)).*(solar'*PV') + CRC(2:T+1,:).*G + GP.*D, 1)  + ...
                RC(1)*PV' + CRC(1,:).*GE' + ITC*C') )
    subject to
        PV >= range1(1);
        PV <= range1(2);
        b1 >= 0;
        b2 >= 0;
        if ds == 0
            b1 + b2 == b_flat;
        end
        D >= 0;
     
        G >= 0;
        G <= ones(T,1)*GE';
        GE >= range2(1);
        GE <= range2(2);
       
        if(N_y == 1)
            C/PP >= max(a_houston_server +  sum(b2,1)'/bu/PP);
            if off == 1
                idle_power = (a_houston_server + sum(b2,1)'/bu/PP)*IP;
            else
                idle_power = C/PP*IP*ones(T,N_y);
            end
            sum(A.*b1,2) + sum(A.*b2,2) == BS;
            sum(b1,2) + sum(b2,2) == BS;
            
            (idle_power+a_houston_power + ...
            sum(b1*(PP-IP)/PP,1)' + ...
                (sum(b2,1)'*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) ...
                    <= D + solar'*PV' + G;     
                
            (idle_power+a_houston_power + ...
                sum(b1*(PP-IP)/PP,1)' + ...
                    (sum(b2,1)'*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) <= C/OP;
            
            if nz == 1
                sum((idle_power+a_houston_power + ...
                    sum(b1*(PP-IP)/PP,1)' + ...
                        (sum(b2,1)'*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) ) <= sum(solar'*PV' + G);
            end
            
            sum(b1,1)' <= a_remaining;
            
            
            D_cap = C/OP;        
            D <= D_cap;   
        else
            C'/PP >= max(a_houston_server +  squeeze(sum(b2,1))/bu/PP);
            if off == 1            
                idle_power = (a_houston_server + squeeze(sum(b2,1))/bu/PP)*IP;   
            else
                idle_power = C/PP*IP*ones(T,N_y);    
            end  
            squeeze(sum(A.*b1,2) + sum(A.*b2,2)) == BS;
            squeeze(sum(b1,2) + sum(b2,2)) == BS;
            
            (idle_power+a_houston_power + ...
            squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ...
                    <= D + solar'*PV' + G;
                
            (idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) <= ones(T,1)*C'/OP;

            if nz == 1
                sum((idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ) <= sum(solar'*PV' + G);
            end            
            squeeze(sum(b1,1)) <= a_remaining;
            
            D_cap = C/OP;        
%             D <= D_cap;              
            D <= ones(T,1)*D_cap';
            
            for y = 2:N_y
                PV(y) >= PV(y-1);
                GE(y) >= GE(y-1);
%                 D_cap(y) >= D_cap(y-1);
                C(y) >= C(y-1);
            end
        end        
cvx_end

if strcmp(cvx_status, 'Solved') == 0
    cvx_status
    save('temp/log.mat');
    error('Increase given capacities');
end
status = cvx_status;
cvx_optval;
% max((idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE)
% mean((idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE)
if(N_y == 1)
    total_demand = 24*365/T*sum((idle_power+a_houston_power + ...
                    sum(b1*(PP-IP)/PP,1)' + ...
                        (sum(b2,1)'*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) ) ;
else
    total_demand =  24*365/T*sum((idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) );
end

total_supply = sum(solar'*PV' + G);
[total_demand; total_supply];
D_cap = max(D);

cost = [RC(1)*PV'; 24*365/T*sum(RC(2)*solar'*PV'); CRC(1,:).*GE'; 24*365/T*sum(CRC(2:T+1,:).*G); 24*365/T*sum(GP.*D)];
capacity = [PV'; GE'; D_cap; C'];
emission = [24*365/T*sum(RE*solar'*PV'); 24*365/T*sum(CRE*G); 24*365/T*CO2_grid*sum(D)];
% ratio = [sum(min(PV*solar,(idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE))/total_demand, sum(G')/total_demand, sum(D')/total_demand];

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
FontSize = 12;
yPlot = 5;
if(N_y == 1)
    demand = (idle_power+a_houston_power + ...
                    sum(b1*(PP-IP)/PP,1)' + ...
                        (sum(b2,1)'*(PP-IP)/PP)).*PUE' ;
else
    demand = (idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y));
end

if plt >= 1
    figure;
    bar(1/24:1/24:T/24,[PV(yPlot)*solar;G(:,yPlot)';D(:,yPlot)']','stacked')
    hold on;
    %plot(1/24:1/24:T/24,(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*PUE,'r',1/24:1/24:T/24,C/OP*ones(1,T),'k', 'LineWidth', 2)
    plot(1/24:1/24:T/24,demand(:,yPlot),'r', 'LineWidth', 2)
    xlabel('month');
    ylabel('power (kW)');
    %legend('PV generation','Gas Engine generation','Grid power','Power demand','Power capacity', 'FontSize', FontSize,'Location', 'northwest');
    legend('PV generation','GE generation','Grid power','Power demand', 'FontSize', FontSize,'Location', 'northwest');
    ylim([0,C/OP*1.1]);
    xlim([0,T/24]);
    set(gca,'XTick',[0.5:1:11.5], 'FontSize', FontSize);
    set(gca,'xticklabel',{'1','2','3','4','5','6','7','8','9','10','11','12'}, 'FontSize', FontSize);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', strcat(location,'.eps'));
%     eps2pdf(strcat(location,'.eps'),'/usr/local/bin/gs');
end

if plt >= 2
    figure;
    power = [idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)+(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*(PUE-ones(1,T))];
    plot(1/24:1/24:T/24, power,'k-','LineWidth',2)
    ylabel('generation (kW)')
    ylim([0,C/OP*1.1])
    set(gca,'XTick',[0.5:1:11.5], 'FontSize', 8);
    set(gca,'xticklabel',{'1','2','3','4','5','6','7','8','9','10','11','12'}, 'FontSize', FontSize);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', strcat(location,'.eps'));
%     eps2pdf(strcat(location,'.eps'),'/usr/local/bin/gs');
    
    figure;
    pie([sum(min(PV*solar,(idle_power+sum(a_houston_power,1) + sum(b1*(PP-IP)/PP,1) + sum(b2,1)*(PP-IP)/PP).*PUE)), sum(G), sum(D)], {'PV','GE','D_cap'})

end

if plt >= 3
    figure;
    bar([a_houston_server;sum(b2,1)/PP/bu]','stacked')
    hold on;
    plot(1:T,C/PP*ones(1,T),'k', 'LineWidth', 2)
    xlabel('month');
    ylabel('server number');
    xlim([1,T]);
    ylim([0,C/PP*1.1]);
end

if plt >= 4
    figure;
    bar([idle_power;sum(a_houston_power,1);sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1);(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*(PUE-ones(1,T))]','stacked')
    hold on;
    plot(1:T,C*ones(1,T),'k',1:T,PV*solar+G','r', 'LineWidth', 2)
    xlabel('month');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C*1.5]);
    legend('idle power','delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable', 'FontSize', FontSize)
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', strcat(location,'.eps'));
%     saveas(gcf,strcat(location,'.fig'))
end

if plt >= 5
    figure;
    bar([(idle_power+sum(a_houston_power,1)+sum(b1/PP*(PP-IP),1)+sum(b2/PP*(PP-IP),1)).*PUE]','stacked')
    hold on;
    plot(1:T,C/OP*ones(1,T),'k',1:T,PV*solar+G','r', 'LineWidth', 2)
    xlabel('month');
    ylabel('kW');
    xlim([1,T]);
    ylim([0,C*1.5]);
    legend('Power demand','Power capaciity','Local generation', 'FontSize', FontSize)
    set(gca,'XTick',[0:6:TS], 'FontSize', 8);
    set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.0 0 4.0 3.5]);
    print ('-depsc', strcat(location,'2.eps'));
%     saveas(gcf,strcat(location,'2.fig'))
end

end

