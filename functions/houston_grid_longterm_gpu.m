function [cost,capacity,emission, total_demand, total_supply] = houston_grid_longterm_gpu(interchange, colocate,nz,ds,off, ...
    T, N_y, ...
    IP,PP,OP, ...
    ITR,ITC, ...
    ITR_gpu,ITC_gpu, alpha, beta,...
    CO2_grid, a, a_gpu, au, ...
    con,BS, BS_gpu,A, A_gpu, S,E,bu,PUE,solar, ...
    range1,RC,RE,SR,range2,CRC,CRE,P,GP,RP,BP,plt,location)

% idle_power = C/PP*IP*ones(T,N_y);
a_houston_server = ceil(a./(au'*ones(T,N_y))/PP); % a_houston_server T,N_y
a_microsoft_server= ceil(a_gpu./(au'*ones(T,N_y))/PP); % a_houston_server T,N_y
a_houston_power = a_houston_server*((PP-IP)*au);
a_microsoft_power = a_microsoft_server*((PP-IP)*au);
if colocate == 1
    a_remaining = a_houston_server*PP*(con-au);
    a_remaining_gpu = a_microsoft_server*PP*(con-au);
else
    a_remaining = zeros(T,N_y);
    a_remaining_gpu = zeros(T,N_y);
end
for y=1:N_y
    b_flat(:,:,y) = BS(:,y)./sum(A(:,:,y),2)*ones(1, T).*A(:,:,y); % todo: repeat A -> (BN, T, N_y); % 
    b_flat_g(:,:,y) = BS_gpu(:,y)./sum(A_gpu(:,:,y),2)*ones(1, T).*A_gpu(:,:,y); % todo: repeat A -> (BN, T, N_y); % 
end

if and(range2(1) > 0 ,range2(1)==range2(2) )
    D_cap = 0;
end
% PV: PV capacity
% GE: Energy storage capacity
% D: Demand capacity
% C: IT capacity
% C_gpu: GPU capacity
cvx_begin
    variables PV(N_y) GE(N_y) G(T,N_y) b1(size(BS,1),T,N_y) b2(size(BS,1),T,N_y) b1_g(size(BS,1),T,N_y) b2_g(size(BS,1),T,N_y) d_g(T,N_y) z1_g(size(BS,1),T,N_y) z2_g(size(BS,1),T,N_y) D(T,N_y) C(N_y) C_gpu(N_y)
    minimize (sum(24*365/T*sum((ones(T,1)*RC(2,:)).*(solar'*PV') + CRC(2:T+1,:).*G + GP.*D, 1)  + ...
                RC(1)*PV' + CRC(1,:).*GE' + ITC*C' + ITC_gpu*C_gpu') )
    subject to
        PV >= range1(1);
%         PV <= range1(2);
        b1 >= 0;
        b2 >= 0;
        b1_g >= 0;
        b2_g >= 0;
        if interchange > 0
          z1_g >= 0;
          z1_g <= b1_g;
          z2_g >= 0;
          z2_g <= b2_g;   
          d_g >= 0;
          d_g <= a_microsoft_server;
        else
          z1_g == 0;          
          z2_g == 0;          
          d_g  == 0;
        end
        
        if ds == 0
            b1 + b2 == b_flat;
            b1_g + b2_g == b_flat;
        end
        D >= 0;
     
        G >= 0;
%         G <= ones(T,1)*GE';
        GE >= range2(1);
%         GE <= range2(2);
        C >= ITR(1);
%         C <= ITR(2);
        C_gpu >= ITR_gpu(1);
%         C_gpu <= ITR_gpu(2);
        
        if(N_y == 1)
            C/PP >= max(a_houston_server +  sum(b2,1)'/bu/PP     + alpha*d_g +  beta*sum(z2_g,1)'/bu/PP);
            C_gpu/PP >= max(a_microsoft_server +  sum(b2_g,1)'/bu/PP - d_g -  sum(z2_g,1)'/bu/PP); % gpu
            if off == 1
                idle_power     = (a_houston_server + sum(b2,1)'/bu/PP + alpha*d_g +  beta*sum(z2_g,1)'/bu/PP )*IP;
                idle_power_gpu = (a_microsoft_server + sum(b2_g,1)'/bu/PP  - d_g -  sum(z2_g,1)'/bu/PP )*IP ; % gpu
            else
                idle_power = (C)/PP*IP*ones(T,N_y); 
                idle_power_gpu = (C_gpu)/PP*IP*ones(T,N_y); % gpu
            end
            if ds ~= 0
              sum(A.*b1,2) + sum(A.*b2,2) == BS;
              sum(b1,2) + sum(b2,2) == BS;

              sum(A.*b1_g,2) + sum(A.*b2_g,2) == BS_gpu; % gpu
              sum(b1_g,2) + sum(b2_g,2) == BS_gpu;
            end
            
            (idle_power+a_houston_power + ...
            sum(b1*(PP-IP)/PP,1)' + ...
                (sum(b2,1)'*(PP-IP)/PP) + ...
            + idle_power_gpu+a_microsoft_power + ... % gpu
            sum(b1_g*(PP-IP)/PP,1)' + ...
                (sum(b2_g,1)'*(PP-IP)/PP)  ...
            + (alpha-1)*d_g + (beta-1)*sum(z1_g*(PP-IP)/PP,1)' + ... % exchanged
                (sum(z2_g,1)'*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) ...    
                    <= D + solar'*PV' + G;     
                
            (idle_power+a_houston_power + ...
                sum(b1*(PP-IP)/PP,1)' + ...
                    (sum(b2,1)'*(PP-IP)/PP) + ...
                     alpha*d_g + beta*(sum(z1_g*(PP-IP)/PP,1)' + beta*(sum(z2_g,1)'*(PP-IP)/PP)) )...
                     .*(PUE'*ones(1,N_y)) <= C/OP;
                  
            (idle_power_gpu+a_microsoft_power + ... % gpu
                sum(b1_g*(PP-IP)/PP,1)' + ...
                    (sum(b2_g,1)'*(PP-IP)/PP)  ...
                    - d_g - (sum(z1_g*(PP-IP)/PP,1)' - (sum(z2_g,1)'*(PP-IP)/PP)) )...
                      .*(PUE'*ones(1,N_y)) <= C_gpu/OP;
            
            if nz == 1
                sum((idle_power+a_houston_power + ...
                    sum(b1*(PP-IP)/PP,1)' + ...
                        (sum(b2,1)'*(PP-IP)/PP)) + ...                      
                (idle_power_gpu+a_microsoft_power + ... % gpu
                    sum(b1_g*(PP-IP)/PP,1)' + ...
                        (sum(b2_g,1)'*(PP-IP)/PP)) + ...
                (alpha-1)*d_g + (beta-1)*(sum(z1_g*(PP-IP)/PP,1)' + (sum(z2_g,1)'*(PP-IP)/PP)) )...                        
                          .*(PUE'*ones(1,N_y)) <= sum(solar'*PV' + G);      
            end
            
            sum(b1,1)' <= a_remaining;            
            sum(b1_g,1)' <= a_remaining_gpu;
            
            D_cap = (C + C_gpu)/OP;        
            D <= D_cap;   
        else
            C'/PP >= max(a_houston_server +  squeeze(sum(b2,1))/bu/PP + alpha*d_g +  beta*squeeze(sum(z2_g,1))/bu/PP );
            C_gpu'/PP >= max(a_microsoft_server +  squeeze(sum(b2_g,1))/bu/PP  - d_g -  squeeze(sum(z2_g,1))/bu/PP ); % gpu            
            if off == 1            
                idle_power = (a_houston_server + squeeze(sum(b2,1))/bu/PP + alpha*d_g +  beta*squeeze(sum(z2_g,1)) ) *IP;   
                idle_power_gpu = (a_microsoft_server + squeeze(sum(b2_g,1))/bu/PP - d_g -  squeeze(sum(z2_g,1)) )*IP ; % gpu
            else
                idle_power = C/PP*IP*ones(T,N_y);    
                idle_power_gpu = C_gpu/PP*IP*ones(T,N_y);  % gpu
            end
            if ds ~= 0
              squeeze(sum(A.*b1,2) + sum(A.*b2,2)) == BS;
              squeeze(sum(b1,2) + sum(b2,2)) == BS;

              squeeze(sum(A.*b1_g,2) + sum(A.*b2_g,2)) == BS_gpu;
              squeeze(sum(b1_g,2)    + sum(b2_g,2))    == BS_gpu;
            end
            
%             (idle_power+a_houston_power + ...
%             squeeze(sum(b1*(PP-IP)/PP,1)) + ...
%                 squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ...
%                     <= D + solar'*PV' + G;
                  
            (idle_power+a_houston_power + ...
            squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                squeeze(sum(b2,1)*(PP-IP)/PP) + ...
            + idle_power_gpu+a_microsoft_power + ... % gpu
            squeeze(sum(b1_g*(PP-IP)/PP,1)) + ...
                squeeze(sum(b2_g,1)*(PP-IP)/PP)  ...
            + (alpha-1)*d_g + (beta-1)*squeeze(sum(z1_g*(PP-IP)/PP,1)) + ... % exchanged
                squeeze(sum(z2_g,1)*(PP-IP)/PP)).*(PUE'*ones(1,N_y)) ...    
                    <= D + solar'*PV' + G;       
                
            (idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP)) ...
                    + alpha*d_g + beta*squeeze(sum(z1_g*(PP-IP)/PP,1)) + beta*squeeze(sum(z2_g,1)*(PP-IP)/PP) ) ...
                    .*(PUE'*ones(1,N_y)) <= ones(T,1)*C'/OP;
                  
            (idle_power_gpu + a_microsoft_power + ... % gpu
                squeeze(sum(b1_g*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2_g,1)*(PP-IP)/PP))  ...
                    - d_g - squeeze(sum(z1_g*(PP-IP)/PP,1)) - squeeze((sum(z2_g,1)*(PP-IP)/PP)) )...
                      .*(PUE'*ones(1,N_y)) <= ones(T,1)*C_gpu'/OP;      

            if nz == 1
                sum((idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ) <= sum(solar'*PV' + G);
                  
                sum((idle_power+a_houston_power + ...
                    squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                        squeeze((sum(b2,1)'*(PP-IP)/PP)) + ...                      
                (idle_power_gpu+a_microsoft_power + ... % gpu
                    squeeze(sum(b1_g*(PP-IP)/PP,1)) + ...
                        squeeze((sum(b2_g,1)'*(PP-IP)/PP))) + ...
                (alpha-1)*d_g + (beta-1)*squeeze(sum(z1_g*(PP-IP)/PP,1) +squeeze(sum(z2_g,1)*(PP-IP)/PP)) )...                        
                          .*(PUE'*ones(1,N_y)) ) <= sum(solar'*PV' + G); 
            end            
            squeeze(sum(b1,1)) <= a_remaining;
            squeeze(sum(b1_g,1)) <= a_remaining_gpu;
            
            D_cap = C/OP;        
%             D <= D_cap;              
            D <= ones(T,1)*D_cap';
            
            for y = 2:N_y
                PV(y) >= PV(y-1);
                GE(y) >= GE(y-1);
%                 D_cap(y) >= D_cap(y-1);
                C(y) >= C(y-1);
                C_gpu(y) >= C_gpu(y-1);
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
D_cap = max(full(D));
% C_gpu
% D_cap
cost = [RC(1)*PV'; 24*365/T*sum(RC(2)*solar'*PV'); CRC(1,:).*GE'; 24*365/T*sum(CRC(2:T+1,:).*G); 24*365/T*sum(GP.*D); ITC*C'; ITC_gpu*C_gpu'];
capacity = [PV'; GE'; D_cap; C'; C_gpu'];
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

