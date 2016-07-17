function [A,BS,S,E] = batch_job_generator(T,BN,O1,A1,B1,O2,A2,B2,BM)
S = ceil(sort(random('Uniform',1,T-1,BN,1)));
BS_raw = random(O2,A2,B2,BN,1);
E = S + ceil(random(O1,A1,B1));
A = zeros(BN, T);
for i = 1:1:BN
%     if E(i) > T
%         BS_raw(i) = BS_raw(i)*(T-S(i))/(E(i)-S(i));
%         E(i) = T;
%     end
%     A(i,S(i):E(i)) = ones(1,E(i)-S(i)+1);
    if E(i) > T
        A(i,S(i):T)  = ones(1,T-S(i)+1);
        A(i,1:E(i)-T) = ones(1,E(i)-T);
    else
        A(i,S(i):E(i)) = ones(1,E(i)-S(i)+1);
    end
end

BS = BS_raw/(mean(BS_raw))*BM*T/BN;