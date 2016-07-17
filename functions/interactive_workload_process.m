function a = interactive_workload_process(location, T, points_per_hour, col)
raw = load(location);
for t = 1:1:T % sum to each hour
    a(t) = sum(raw((t-1)*points_per_hour+1:t*points_per_hour,col));
end