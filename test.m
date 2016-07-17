load('temp/log.mat');
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
    variables PV(N_y) GE(N_y) G(T,N_y) b1(size(BS,1),T,N_y) b2(size(BS,1),T,N_y) D(T,N_y)
    minimize (sum(24*365/T*sum((ones(T,1)*RC(2,:)).*(solar'*PV') + CRC(2:T+1,:).*G + GP.*D, 1)  + ...
                RC(1)*PV' + CRC(1,:).*GE'))
    subject to
        PV >= range1(1);
        PV <= range1(2);
        b1 >= 0;
        b2 >= 0;
        if ds == 0
            b1 + b2 == b_flat;
        end
        D >= 0;
        D <= D_cap;        
        G >= 0;
        G <= ones(T,1)*GE';
        GE >= range2(1);
        GE <= range2(2);
       
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
            
%             D <= D_cap;
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
            
%             D <= ones(T,1)*D_cap';
            for y = 2:N_y
                PV(y) >= PV(y-1);
                GE(y) >= GE(y-1);
%                 D_cap(y) >= D_cap(y-1);
            end
        end        
cvx_end