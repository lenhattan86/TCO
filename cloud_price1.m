clear;
clc;

% total capacity, cooling
C = 1000;
PUE = 1.3;
T = 24;

% interactive workload
%a = interactive_process('traces\SAPnew\sapTrace.tab', T, 12, 4, 2, 0.5*C);
a = zeros(1,T);

% batch jobs
%bid_high = 12;
%bid_low = 3;
bid_high = 10;
bid_low = 4;
step_size = 0.1;
batch_size = 0.4;
batch_arrival = batch_size*ceil(random('Uniform',0,20,(bid_high-bid_low)/step_size,T));

% PV generation
R = trace_process('traces\solar-one-week.csv', T, 12, 1, 3, C*PUE);
%R = zeros(1,T);

% electricity price
p_raw = load('traces\pa-electric-price.txt'); % electricity price trace
PS = 0.1; % scale to cents/kWh
ep = p_raw*PS; % scale to cents/kWh
%ep = mean(ep)*ones(1,T);
oc = 1000;

% solver1
batch_remaining = zeros(size(batch_arrival,1),T);
for t = 1:1:T
    if t == 1
        batch_curr(:,t) = batch_arrival(:,t);
    else
        batch_curr(:,t) = batch_remaining(:,t-1) + batch_arrival(:,t);
    end
    profit_opt = -oc;
    price_opt = bid_low;
    revenue_opt = 0;
    cost_opt = 0;
    demand_opt = 0;
    for p = bid_high-step_size:-step_size:bid_low
        demand_p = sum(batch_curr(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
        if demand_p + a(t) > C
            break;
        else
            profit_p = demand_p*p - ep(t)*max(0,demand_p*PUE-max(0,R(t)-a(t)*PUE))-oc;
            cost_p = ep(t)*max(0,demand_p*PUE-max(0,R(t)-a(t)*PUE));
            revenue_p = demand_p*p;
            if profit_p > profit_opt
                profit_opt = profit_p;
                cost_opt = cost_p;
                revenue_opt = revenue_p;
                demand_opt = demand_p;
                price_opt = p;
            end                
        end
    end
    price(t) = price_opt;
    profit(t) = profit_opt;
    revenue(t) = revenue_opt;
    cost(t) = cost_opt;
    demand(t) = demand_opt;
    batch_finish(:,t) = [zeros(floor((price_opt+step_size/10-bid_low)/step_size),1);batch_curr(ceil((price_opt+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining(:,t) = batch_curr(:,t) - batch_finish(:,t);
end

total_profit = sum(profit);
total_revenue = sum(revenue);
total_cost = sum(cost);
total_oc = oc*T;
total_demand = sum(demand);
avg_price = mean(price);

figure;
bar([a',sum(batch_finish,1)',(a'+sum(batch_finish,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1.eps');
saveas(gcf,'results/cloud1.fig')

figure;
bar(batch_finish','stacked')

%{
figure;
plot(1:T,price)
xlim([1,T]);
xlabel('hour');
ylabel('price (cents/kWh)');
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/price1.eps');
saveas(gcf,'results/price1.fig')

figure;
plot(1:T,profit);

%}


% solver2: optimal flat price
p_opt_flat = bid_high;
total_profit_opt_flat = 0;
for p = bid_high-step_size:-step_size:bid_low
    batch_remaining_flat = zeros(size(batch_arrival,1),T);
    feasible = 1;
    for t = 1:1:T
        if t == 1
            batch_curr_flat(:,t) = batch_arrival(:,t);
        else
            batch_curr_flat(:,t) = batch_remaining_flat(:,t-1) + batch_arrival(:,t);
        end
        demand_p_flat = sum(batch_curr_flat(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
        if demand_p_flat + a(t) > C
            feasible = 0;
            break;
        else
            profit_flat(t) = demand_p_flat*p - ep(t)*max(0,demand_p_flat*PUE-max(0,R(t)-a(t)*PUE))-oc;              
        end
    end
    if feasible == 1
        total_profit_flat = sum(profit_flat);
        if (total_profit_flat > total_profit_opt_flat)
             p_opt_flat = p
             total_profit_opt_flat = total_profit_flat
        end
    end
end

batch_remaining_flat = zeros(size(batch_arrival,1),T);
p = p_opt_flat;
for t = 1:1:T
    if t == 1
        batch_curr_flat(:,t) = batch_arrival(:,t);
    else
        batch_curr_flat(:,t) = batch_remaining_flat(:,t-1) + batch_arrival(:,t);
    end
    demand_p_flat = sum(batch_curr_flat(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
    profit_flat(t) = demand_p_flat*p - ep(t)*max(0,demand_p_flat*PUE-max(0,R(t)-a(t)*PUE))-oc;
    revenue_flat(t) = demand_p_flat*p;
    cost_flat(t) = ep(t)*max(0,demand_p_flat*PUE-max(0,R(t)-a(t)*PUE));
    demand_flat(t) = demand_p_flat;
    batch_finish_flat(:,t) = [zeros(floor((p+step_size/10-bid_low)/step_size),1);batch_curr_flat(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining_flat(:,t) = batch_curr_flat(:,t) - batch_finish_flat(:,t);
end
total_profit_flat = sum(profit_flat);
total_revenue_flat = sum(revenue_flat);
total_cost_flat = sum(cost_flat);
total_oc_flat = oc*T;
total_demand_flat = sum(demand_flat);
avg_price_flat = mean(p);
price_flat = p*ones(1,T);

figure;
bar([a',sum(batch_finish_flat,1)',(a'+sum(batch_finish_flat,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_flat.eps');
saveas(gcf,'results/cloud1_flat.fig')

figure;
bar(batch_finish_flat','stacked')

%{
figure;
bar([a',sum(batch_finish_flat,1)',(a'+sum(batch_finish_flat,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1.eps');
saveas(gcf,'results/cloud1.fig')



figure;
plot(1:T,profit,1:T,profit_flat);

figure;
bar(batch_finish_flat','stacked')


figure;
bar([total_profit,total_profit_flat]);

%}

% solver3: low flat price
batch_remaining_low = zeros(size(batch_arrival,1),T);
p = max(p_opt_flat-step_size*10,bid_low);
for t = 1:1:T
    if t == 1
        batch_curr_low(:,t) = batch_arrival(:,t);
    else
        batch_curr_low(:,t) = batch_remaining_low(:,t-1) + batch_arrival(:,t);
    end
    demand_p_low = sum(batch_curr_low(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
    profit_low(t) = demand_p_low*p - ep(t)*max(0,demand_p_low*PUE-max(0,R(t)-a(t)*PUE))-oc;
    revenue_low(t) = demand_p_low*p;
    cost_low(t) = ep(t)*max(0,demand_p_low*PUE-max(0,R(t)-a(t)*PUE));
    demand_low(t) = demand_p_low;
    batch_finish_low(:,t) = [zeros(floor((p+step_size/10-bid_low)/step_size),1);batch_curr_low(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining_low(:,t) = batch_curr_low(:,t) - batch_finish_low(:,t);
end
total_profit_low = sum(profit_low);
total_revenue_low = sum(revenue_low);
total_cost_low = sum(cost_low);
total_oc_low = oc*T;
total_demand_low = sum(demand_low);
avg_price_low = mean(p);
price_low = p*ones(1,T);

figure;
bar([a',sum(batch_finish_low,1)',(a'+sum(batch_finish_low,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_low.eps');
saveas(gcf,'results/cloud1_low.fig')

figure;
bar(batch_finish_low','stacked')

% solver4: high flat price
batch_remaining_high = zeros(size(batch_arrival,1),T);
p = min(p_opt_flat+step_size*10,bid_high);
for t = 1:1:T
    if t == 1
        batch_curr_high(:,t) = batch_arrival(:,t);
    else
        batch_curr_high(:,t) = batch_remaining_high(:,t-1) + batch_arrival(:,t);
    end
    demand_p_high = sum(batch_curr_high(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
    profit_high(t) = demand_p_high*p - ep(t)*max(0,demand_p_high*PUE-max(0,R(t)-a(t)*PUE))-oc;
    revenue_high(t) = demand_p_high*p;
    cost_high(t) = ep(t)*max(0,demand_p_high*PUE-max(0,R(t)-a(t)*PUE));
    demand_high(t) = demand_p_high;
    batch_finish_high(:,t) = [zeros(floor((p+step_size/10-bid_low)/step_size),1);batch_curr_high(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining_high(:,t) = batch_curr_high(:,t) - batch_finish_high(:,t);
end
total_profit_high = sum(profit_high);
total_revenue_high = sum(revenue_high);
total_cost_high = sum(cost_high);
total_oc_high = oc*T;
total_demand_high = sum(demand_high);
avg_price_high = mean(p);
price_high = p*ones(1,T);

figure;
bar([a',sum(batch_finish_high,1)',(a'+sum(batch_finish_high,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_high.eps');
saveas(gcf,'results/cloud1_high.fig')

figure;
bar(batch_finish_high','stacked')


% solver4: binary price
batch_remaining_binary = zeros(size(batch_arrival,1),T);
p(1:8) = mean(price(1:8));
p(9:16) = mean(price(9:16));
p(17:24) = mean(price(17:24));

for t = 1:1:T
    if t == 1
        batch_curr_binary(:,t) = batch_arrival(:,t);
    else
        batch_curr_binary(:,t) = batch_remaining_binary(:,t-1) + batch_arrival(:,t);
    end
    demand_p_binary = sum(batch_curr_binary(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
    profit_binary(t) = demand_p_binary*p(t) - ep(t)*max(0,demand_p_binary*PUE-max(0,R(t)-a(t)*PUE))-oc;
    revenue_binary(t) = demand_p_binary*p(t);
    cost_binary(t) = ep(t)*max(0,demand_p_binary*PUE-max(0,R(t)-a(t)*PUE));
    demand_binary(t) = demand_p_binary;
    batch_finish_binary(:,t) = [zeros(floor((p(t)+step_size/10-bid_low)/step_size),1);batch_curr_binary(ceil((p(t)+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining_binary(:,t) = batch_curr_binary(:,t) - batch_finish_binary(:,t);
end
total_profit_binary = sum(profit_binary);
total_revenue_binary = sum(revenue_binary);
total_cost_binary = sum(cost_binary);
total_oc_binary = oc*T;
total_demand_binary = sum(demand_binary);
avg_price_binary = mean(p);
price_binary = p;

figure;
bar([a',sum(batch_finish_binary,1)',(a'+sum(batch_finish_binary,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_binary.eps');
saveas(gcf,'results/cloud1_binary.fig')

figure;
bar(batch_finish_binary','stacked')

% solver5: maximize utilization
batch_remaining_mu = zeros(size(batch_arrival,1),T);
for t = 1:1:T
    if t == 1
        batch_curr_mu(:,t) = batch_arrival(:,t);
    else
        batch_curr_mu(:,t) = batch_remaining_mu(:,t-1) + batch_arrival(:,t);
    end
    profit_opt_mu = -oc;
    price_opt_mu = bid_low;
    revenue_opt_mu = 0;
    cost_opt_mu = 0;
    demand_opt_mu = 0;
    for p = bid_high-step_size:-step_size:bid_low
        demand_p_mu = sum(batch_curr_mu(ceil((p+step_size/10-bid_low)/step_size):size(batch_arrival,1),t));
        if demand_p_mu + a(t) > C || p < ep(t)
            break;
        else
            profit_p_mu = demand_p_mu*p - ep(t)*max(0,demand_p_mu*PUE-max(0,R(t)-a(t)*PUE))-oc;
            cost_p_mu = ep(t)*max(0,demand_p_mu*PUE-max(0,R(t)-a(t)*PUE));
            revenue_p_mu = demand_p_mu*p;
            if demand_p_mu > demand_opt_mu
                profit_opt_mu = profit_p_mu;
                cost_opt_mu = cost_p_mu;
                revenue_opt_mu = revenue_p_mu;
                demand_opt_mu = demand_p_mu;
                price_opt_mu = p;
            end                
        end
    end
    price_mu(t) = price_opt_mu;
    profit_mu(t) = profit_opt_mu;
    revenue_mu(t) = revenue_opt_mu;
    cost_mu(t) = cost_opt_mu;
    demand_mu(t) = demand_opt_mu;
    batch_finish_mu(:,t) = [zeros(floor((price_opt_mu+step_size/10-bid_low)/step_size),1);batch_curr_mu(ceil((price_opt_mu+step_size/10-bid_low)/step_size):size(batch_arrival,1),t)]';
    batch_remaining_mu(:,t) = batch_curr_mu(:,t) - batch_finish_mu(:,t);
end

total_profit_mu = sum(profit_mu);
total_revenue_mu = sum(revenue_mu);
total_cost_mu = sum(cost_mu);
total_oc_mu = oc*T;
total_demand_mu = sum(demand_mu);
avg_price_mu = mean(price_mu);

figure;
bar([a',sum(batch_finish_mu,1)',(a'+sum(batch_finish_mu,1)')*(PUE-1)],'stacked');
hold on;
plot(1:T,C*ones(1,T),'k',1:T,R,'r', 'LineWidth', 2);
xlim([1,T]);
xlabel('hour');
ylabel('kW');
legend('delay-sensitive workload','delay-tolerant workload','cooling power','IT capacity','renewable')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_mu.eps');
saveas(gcf,'results/cloud1_mu.fig')

figure;
bar(batch_finish_mu','stacked')


figure;
plot(1:T, (a'+sum(batch_finish,1)')*PUE, 'k-', 1:T, (a'+sum(batch_finish_flat,1)')*PUE, 'b-', 1:T, (a'+sum(batch_finish_low,1)')*PUE, 'b--', 1:T, (a'+sum(batch_finish_high,1)')*PUE, 'b-.',...
    1:T, (a'+sum(batch_finish_binary,1)')*PUE, 'k--', 1:T, (a'+sum(batch_finish_mu,1)')*PUE, 'r-')
xlim([1,T]);
ylim([0,C*PUE*1.2]);
xlabel('hour');
ylabel('kW');
legend('dynamic price','optimal flat price','low flat price','high flat price','binary flat price','maximize utilization')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_demand.eps');
saveas(gcf,'results/cloud1_demand.fig')

figure;
plot(1:T, price, 'k-', 1:T, price_flat, 'b-', 1:T, price_low, 'b--', 1:T, price_high, 'b-.',...
    1:T, price_binary, 'k--', 1:T, price_mu, 'r-')
xlim([1,T]);
ylim([bid_low,bid_high]);
xlabel('hour');
ylabel('price (cents/kWh)');
legend('dynamic price','optimal flat price','low flat price','high flat price','binary flat price','maximize utilization')
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 5.0 2.8]);
print ('-depsc', 'results/cloud1_price.eps');
saveas(gcf,'results/cloud1_price.fig')

csvwrite('results\price1_summary.csv',[total_revenue,total_cost,total_oc,total_profit,total_demand,avg_price;...
    total_revenue_flat,total_cost_flat,total_oc_flat,total_profit_flat,total_demand_flat,avg_price_flat;...
    total_revenue_low,total_cost_low,total_oc_low,total_profit_low,total_demand_low,avg_price_low;...
    total_revenue_high,total_cost_high,total_oc_high,total_profit_high,total_demand_high,avg_price_high;...
    total_revenue_binary,total_cost_binary,total_oc_binary,total_profit_binary,total_demand_binary,avg_price_binary;...
    total_revenue_mu,total_cost_mu,total_oc_mu,total_profit_mu,total_demand_mu,avg_price_mu]);


% future: gas engine, co-location, batch job size, optimal