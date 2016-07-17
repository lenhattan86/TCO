
R_raw = zeros(12,24*31);
PVN = zeros(12,1);
CF = zeros(12,1);
RCF = load('traces\PV_CF.csv');
for i = 1:1:12
    tmp = load(strcat('traces/Houston/',int2str(i),'.csv'));
    PVN(i) = size(tmp,1)/12;
    for t = 1:1:PVN(i)
        R_raw(i,t) = max(0,mean(tmp((t-1)*12+1:(t-1)*12+12,6)));
    end
end

for i = 1:1:12
    CF(i) = mean(R_raw(i,:))/max(max(R_raw));
end

for i = 1:1:12
    R((i-1)*24+1:(i-1)*24+24) = sum(reshape(R_raw(i,1:PVN(i)),24,PVN(i)/24),2);
end

R = R/max(R);