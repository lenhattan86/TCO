function [cost, emission, total_demand, total_supply] = realize_houston_grid_longterm(colocate,nz,ds,off, ...
    T, N_y, ...
    IP,PP,OP,C,CO2_grid,a,au, ...
    con,BS,A,S,E,bu,PUE,solar, ...
    range1,RC,RE,SR,range2,CRC,CRE,P,GP,RP,BP, capacity)

PV = capacity(1,:)';
GE = capacity(2,:)';

idle_power = C/PP*IP*ones(T,N_y);
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
D_cap = C/OP;
if and(range2(1) > 0 ,range2(1)==range2(2) )
    D_cap = 0;
end
% PV: PV capacity
% GE: Energy storage capacity
% D: Demand capacity
cvx_begin
    variables G(T,N_y) b1(size(BS,1),T,N_y) b2(size(BS,1),T,N_y) D(T,N_y)
    minimize (sum(24*365/T*sum((ones(T,1)*RC(2,:)).*(solar'*PV') + CRC(2:T+1,:).*G + GP.*D, 1)  + ...
                RC(1)*PV' + CRC(1,:).*GE'))
    subject to
        b1 >= 0;
        b2 >= 0;
        if ds == 0
            b1 + b2 == b_flat;
        end
        D >= 0;
        D <= D_cap;        
        G >= 0;
        G <= ones(T,1)*GE';
       
        if(N_y == 1)
            C/PP >= max(a_houston_server +  sum(b2,1)'/bu/PP);
            if off == 1
                idle_power = (a_houston_server + sum(b2,1)'/bu/PP)*IP;
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

        else
            C/PP >= max(a_houston_server +  squeeze(sum(b2,1))/bu/PP);
            if off == 1            
                idle_power = (a_houston_server + squeeze(sum(b2,1))/bu/PP)*IP;            
            end  
            squeeze(sum(A.*b1,2) + sum(A.*b2,2)) == BS;
            squeeze(sum(b1,2) + sum(b2,2)) == BS;
            
            (idle_power+a_houston_power + ...
            squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ...
                    <= D + solar'*PV' + G;
                
            (idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) <= C/OP;

            if nz == 1
                sum((idle_power+a_houston_power + ...
                squeeze(sum(b1*(PP-IP)/PP,1)) + ...
                    squeeze((sum(b2,1)*(PP-IP)/PP))).*(PUE'*ones(1,N_y)) ) <= sum(solar'*PV' + G);
            end            
            squeeze(sum(b1,1)) <= a_remaining;
            
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
emission = [24*365/T*sum(RE*solar'*PV'); 24*365/T*sum(CRE*G); 24*365/T*CO2_grid*sum(D)];
end

