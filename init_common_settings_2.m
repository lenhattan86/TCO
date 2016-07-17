addpath('functions/');

has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_solver Gurobi;
cvx_precision low;
cvx_clear;

%fig_path = '../figs/';
fig_path = 'results/';
% general parameters
T = 24*7; % planning period length, in hours
IP = 0.04; % idle power consumption of a server, in kW
PP = 0.2; % peak power consumption of a server, in kW
CO2_grid = 0.5*1000; % g CO2/kWh of grid power consumption
DC_power = 1000; % Data center power capacity, in kW
OP = 5/6; % 0.8 % over-provisioning ratio of IT and power
IT = DC_power*OP;% Data center IT capacity, in kW
con = 0.9; % maximum utilization after consolidation

% interactive workloads
IN = 1; % number of interactive workloads
IL = {'traces\SAPnew\sapTrace.tab'};
a = zeros(IN,T);
au = 0.4; % average utilization of interactive workloads
IPU = 0.7; % interactive workload peak utilization
IM = IT*au*IPU; % peak of interactive workload consumption, in servers
PMR = 2.5; % peak to mean ratio, should be smaller than current PMR
for i = 1:1:IN
    a(i,:) = interactive_process(char(IL(i)), T, 12, 4, PMR, IM(i));
end

% batch job
BN = 1; % average number of batch job arrivals per timeslot
BM = 0.25; % batch job ratio, compared with interactive workload
[A,BS,S,E] = batch_job_generator(T,BN*T,'Uniform',23.99,23.99,'Uniform',1,1,BM/(1-BM)*sum(mean(a,2)./au)/con); % generate batch jobs based on statistical properties
bu = 1; % average utilization of batch job when running alone.

% uncontrolled renewable generation, e.g., PV, wind
RN = 1; % number of renewable generation
RL = {'traces\solar-one-week-07012012.csv'};
R = zeros(IN,T);
RP = 0.3; % capacity factor
RR = [0, 1.0*DC_power]; % renewable capacity range, in kW
for i = 1:1:RN
    R(i,:) = trace_process(char(RL(i)), T, 12, 1, 1, RP(i));
end
RC = [2600/20, 0.005+0]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
RE = 32; % emmision, in gram CO2 eq

% controlled renewable generation, e.g., gas engine
CRN = 1; % number of controlled renewable generation
SR = 1/6; % storage requirement for 1kW, in kWh
CRR = [0,1400]; % controled renewable capacity range, in kW
% CRC = [1000/20+SR*200/4, 0.005+0.06]; % [install cost ($/kW/year), maintenance cost + operational cost ($/kWh)]
CRC = [1000/20+SR*200/4, 0.005+0.06*ones(1,T)]; % Tan added this line to fix compiler error.
CRE = 443; % emmision, in gram CO2 eq

% cooling efficiency
PUE = 1.3*ones(1,T);

% price setting [power, cooling, PV, IT]
P = [600, 800, 3200/20, 1000];
GP = 0.08*ones(1,T); % grid electricity price
BP = 0.04*ones(1,T); % sell back price
RP = 0.14*ones(1,T); % renewable electricity price

% energy storage loss
Loss = [0.99, 1];