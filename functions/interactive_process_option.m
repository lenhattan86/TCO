function a = interactive_process(option,location, T, points_per_hour, col, PMR, var)
raw = load(location);
for t = 1:1:T % sum to each hour
    a_raw(t) = max(0,sum(raw((t-1)*points_per_hour+1:t*points_per_hour,col)));
end
PMR_raw = max(a_raw)/mean(a_raw);
assert(PMR_raw>PMR);
mean_raw = mean(a_raw);
a_raw = a_raw + mean_raw*((PMR_raw-1)/(PMR-1)-1);   
if option == 1
	a = a_raw/max(a_raw)*var;
elseif option == 2
    a = a_raw/mean(a_raw)*var;
end
assert(abs(max(a)/mean(a)-PMR) <= 0.01*PMR)