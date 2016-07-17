function a = trace_process(location, T, points_per_hour, col, option, var)
raw = load(location);
for t = 1:1:T % sum to each hour
    a_raw(t) = max(0,sum(raw((t-1)*points_per_hour+1:t*points_per_hour,col)));
end
if option == 1
    a = a_raw/mean(a_raw)*var;
elseif option == 2
    a = a_raw/var;
elseif option == 3
    a = a_raw/max(a_raw)*var;
elseif option == 4
    a = mean(reshape(a_raw,var,T/var),2)';
end