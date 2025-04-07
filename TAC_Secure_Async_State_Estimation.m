close all; clear; clc
%% This is the secure state estimation simulation on the IEEE 14-bus system
N = 14; % number of buses
n = N*2; % number of states
m = N*3; % number of sensors
Ts = 1e-2; % sampling time
ns = 10; % total time steps
%% Network parameters
Adj = zeros(14,14);
Adj(1,2) = 1;  Adj(1,5) = 1;
Adj(2,3) = 1;  Adj(2,4) = 1;  Adj(2,5) = 1;
Adj(3,4) = 1;
Adj(4,5) = 1;  Adj(4,7) = 1;  Adj(4,9) = 1;
Adj(5,6) = 1;
Adj(6,11)= 1; Adj(6,12) = 1; Adj(6,13) = 1;
Adj(7,8) = 1; Adj(7,9) = 1;
Adj(9,10)= 1; Adj(9,14) = 1; 
Adj(10,11) = 1;
Adj(12,13) = 1;
Adj(13,14) = 1;
Adj = Adj*15;
Adj = Adj + Adj';
M = 100*eye(14);
M(1,1) = 10; M(2,2) = 10; 
M(3,3) = 10;
M(6,6) = 10; 
M(8,8) = 10;
L = diag(sum(Adj,1)) - Adj; % Laplacian matrix representing the network
D = 2*eye(N); % Damping matrix
%% System parameters
A = [zeros(N,N), eye(N);
     -M\L, -M\D]; % Matrix A represents the system dynamics
C0 = [0, 1;
      1, 0;
      0, 1]; 
C = []; % Measurement matrix C
for i = 1:N
    ei = zeros(1,N); ei(i) = 1;
    C = [C;kron(C0,ei)];
end
Q = 5e-3*diag(1+0.1*rand(n,1)); % psd matrix Q for process disturbance
R = 1e-3*diag(1+0.1*rand(m,1)); % psd matrix R for measurement noise
rbar = max(max(R)); % upper bound of covariance Ri
x0 = randn(n,1); % initial state x(0)
Pk = eye(n); % Sigma = Pk(0)

%% Initialization of matrices G, H, W
Gk = cell(m,ns); 
H = cell(m,1);
Hk = [];
ssc = ss(A,zeros(n,1),C,zeros(m,1));
ssd = c2d(ssc,Ts*5,'zoh');
Ad = ssd.A;
for i = 1:m
    Gk{i,1} = zeros(n,n);
    H{i,1} = zeros(n,n);
    Ci = C(i,:);
    Obi = obsv(Ad,Ci);
    r = rank(Obi);
    [~,pvi] = rref(Obi);
    for j = 1:length(pvi)
        Gk{i,1}(pvi(j),pvi(j)) = 1;
        H{i,1}(pvi(j),pvi(j)) = 1; % Matrix H(0)
    end
    Hk = [Hk;H{i,1}];
end
Gsum = zeros(n,n);
for i = 1:m
    Gsum = Gsum + Gk{i,1};
end
iGsum = Gsum\eye(n);
for i = 1:m
    Gk{i,1} = Gk{i,1} * iGsum; % Matrix G(0)
end
Wkb = cell(1,ns);
D = zeros(n*m,n*m); 
for i = 1:m
    for j = 1:m
        if i == j
            D(1+(i-1)*n:i*n,1+(j-1)*n:j*n) = Gk{i,1} + (m-1) * eye(n);
        else
            D(1+(i-1)*n:i*n,1+(j-1)*n:j*n) = - eye(n);
        end
    end
end
Wkb{1,1} = D .* (kron(ones(m,m),Pk)); % Matrix W(0)

%% Variables
x = zeros(n,ns); % true state
x(:,1) = x0; % inital state
xh = zeros(n,ns); % state estimate by KF
xhc = zeros(n,ns); % state estimate oracle by KF
rmse_xh = zeros(1,ns); % relative mean square error of x hat
rmse_xhc = zeros(1,ns); % relative mean square error of x hat oracle
xls = zeros(n,ns); % state by least square fusion 
muls = zeros(n*m,ns); % noise estimation
varls = ones(n*m,ns); % attack estimation
rmse_xls = zeros(1,ns); % relative mean square error of x least square
w = sqrt(Q) * randn(n,ns); % process noise
v = sqrt(R) * randn(m,ns); % measurement noise
zeta = cell(m,1); % all local estimator
gammaLS = 2; % gamma LASSO
for s = 1:m
    zeta{s} = zeros(n,ns);
    zeta{s}(:,1) = Gk{s,1} * x0; % initialize zeta(0)
end
zetasum = zeros(n,ns); % sum of all local estimator
y = zeros(m,ns); % output measurement
yc = zeros(m,ns); % oracle output measurement
phi = zeros(m,ns); % measurement available index
avail_set = [zeros(1,5) 1]; % a list of 0 and 1, this is used to determine whether a sensor sends measurement

%% Asynschronous time variable
tk = 0; % starting samping time
ts = tk; % collection of sampling times
dtk = 0; % time interval between 2 measurements sent

% Attack scenario
flag_fdi = 1; % flag for false data injection attacks
flag_ts = 1; % flag for time-stamp attacks
flag_dos = 1; % flag for denial-of-service attacks
flag_fdg = 1; % flag for false data generation attacks
sc_fdi = 3 * 3 - randi([0 2],1,1); % false data injection attacks on a random sensor of bus 3
sc_ts = 3 * 5   - randi([0 2],1,1); % time-stamp attacks on a random sensor of bus 5
sc_dos = 3 * 4  - randi([0 2],1,1); % denial-of-service attacks on a random sensor of bus 4
sc_fdg = 3 * 2  - randi([0 2],1,1); % fake data generation attacks on a random sensor of bus 2
att_activ = cell(4,1); % collection of times attacks happend

% Attack duration
fdi_start = 5; fdi_end = round(ns); % false data injection attacks start after 5 measurments
ts_start = 6; ts_end = round(ns); % time-stamp attacks start after 6 measurments
dos_start = 7; dos_end = round(ns); % denial-of-service attacks after 7 measurments
fdg_start = 8; fdg_end = round(ns); % fake data generation attacks after 8 measurments

ts_record = []; % recorded time when time-stamp attacks happend

% Optimizer
xlsv = sdpvar(n,1); % optimization variable x for the least square
mu = sdpvar(n*m,1); % optimization variable mu for the least square
var = sdpvar(n*m,1); % optimization variable vartheta for the least square
options = sdpsettings('solver','mosek','verbose',0,'debug',1);

%% Simulation
for k = 2:ns
    % time asynchronous
    dtk = Ts * randi([1 5],1,1); % time interval between 2 consecutive sampling measurements
    tk = tk + dtk; % next samping time
    ts = [ts,tk]; % collection of sampling times
    Ak = expm(A * dtk); % system matrix Ak
    iAk = expm(-A * dtk); % its inverse Ak^{-1}
    Qk = Q * dtk; % process noise covariance matrix Qk
    Rdk = R; % measurement noise covariance matrix Rk 

    % System dynamics
    x(:,k) = Ak * x(:,k-1) + w(:,k) * dtk; % system dynamics
    y(:,k) = C * x(:,k) + v(:,k); % measurement if sensor samples
    yc(:,k) = y(:,k); % oracle measurement without attacks

    % Measurement availability
    for s = 1:m
        phi(s,k) = avail_set(randi([1 length(avail_set)],1,1)); % some sensors sample if phi(k) = 1
    end
    if max(phi(:,k)) < 1
        phi(randi([1 m],1,1),k) = 1; % make sure, at least one sensor samples
    end
    phi(sc_fdi,fdi_start) = 1; % make sure, false data injection attacks happend
    phi(sc_ts,ts_start) = 1; % make sure, time-stamp attacks happend
    phi(sc_dos,dos_start) = 1; % make sure, denial-of-service attacks happend
    phi(sc_fdg,fdg_start) = 1; % make sure, false data generation attacks happend

    if k <= 4
        phi(:,k) = 1; % assume that all the sensors sample measurements at the first 3 steps
    end

    % Adversarial manipulation
    % false data injection attacks
    if flag_fdi == 1 && phi(sc_fdi,k) == 1 && k >= fdi_start && k <= fdi_end 
        y(sc_fdi,k) = y(sc_fdi,k) + 2 + randn(1,1); % inject a random signal to measurement
        att_activ{1,1} = [att_activ{1,1},tk]; % collect times the attack happends
    end

    % time-stamp attacks
    if flag_ts == 1 && k >= ts_start && k <= ts_end
        if phi(sc_ts,k) == 1
            ts_record = [ts_record;y(sc_ts,k)]; % record the measurement
            phi(sc_ts,k) = 0; % disable the measurement sent to fusion center
            att_activ{2,1} = [att_activ{2,1},tk]; % collect times the attack happends
        elseif phi(sc_ts,k) == 0 && ~isempty(ts_record)
            ts_take = randi(length(ts_record)); % take a random true measurement
            y(sc_ts,k) = ts_record(ts_take); % release the true measurment at different time-stamp          
            ts_record(ts_take) = []; % empty the record to make sure that no measurement is duplicated
            phi(sc_ts,k) = 1; % enable the measurement sent to fusion center
        end
    end
    
    % denial-of-service attacks
    if flag_dos == 1 && k >= dos_start && k <= dos_end
        if phi(sc_dos,k) == 1
            phi(sc_dos,k) = 0; % disable the measurement sent to fusion center
            att_activ{3,1} = [att_activ{3,1},tk]; % collect times the attack happends
        end
    end
    
    % fake data generation attacks
    if flag_fdg == 1 && k >= fdg_start && k <= fdg_end
        if phi(sc_fdg,k) == 0
            phi(sc_fdg,k) = avail_set(randi([1 length(avail_set)],1,1)); % enable the measurement sent to fusion center       
        end
        if phi(sc_fdg,k) == 1
            y(sc_fdg,k) = 0.3 + randn(1,1); % generate completely new false data measurements
            att_activ{4,1} = [att_activ{4,1},tk]; % collect times the attack happends
        end
    end

    % State estimation matrices
    Ck = diag(phi(:,k)) * C; % compute matrix Ck by sampling time
    Rk = diag(phi(:,k)) * Rdk; % compute matrix Rk by sampling time
    % Prediction step
    yk = diag(phi(:,k)) * y(:,k);
    xh(:,k) = Ak * xh(:,k-1);
    xhc(:,k) = Ak * xhc(:,k-1);
    Pk = Ak * Pk * Ak' + Qk;
    % Update step
    Kk = Pk * Ck' * pinv(Ck*Pk*Ck' + Rk);
    Pk = (eye(n) - Kk*Ck)*Pk;
    xh(:,k) = xh(:,k) + Kk * (yk - Ck*xh(:,k)); % state estimate by KF with possibly attacked measurements
    xhc(:,k) = xhc(:,k) + Kk * (diag(phi(:,k)) * yc(:,k) - Ck*xhc(:,k)); % state estimate by KF without attacks
    % Relative root mean square
    rmse_xh(k) = norm((xh(:,k) - x(:,k))); % root mean square error of the state estimate with the true state by KF with possibly attacked measurements
    rmse_xhc(k) = norm((xhc(:,k) - x(:,k))); % root mean square error of the state estimate with the true state by KF without attacks

    %% Asynschronous Sensor Fusion
    Pik = Ak - Kk * Ck * Ak; % recursively compute matrix Pi
    % Local estimator
    for s = 1:m
        zeta{s}(:,k) = Pik * zeta{s}(:,k-1) + Kk(:,s) * yk(s); % dyanmics of local estimator zeta_i
        zetasum(:,k) = zetasum(:,k) + zeta{s}(:,k); % compute sum of all local estimator to double-check the state estimate by summation of all the local estimator
    end
    % Least square fusion
    Qkv = [];
    Kkv = [];
    zetav = [];
    Gkv = [];
    Gcheck = zeros(n,n);
    for s = 1:m    
        Gk{s,k} = Pik * Gk{s,k-1} * iAk + Kk(:,s) * Ck(s,:); % dynamics of G_i
        zetav = [zetav;zeta{s}(:,k)]; % stack up the local estimtor as a vector
        Gkv = [Gkv;Gk{s,k}];
        Qkv = [Qkv;Pik*Gk{s,k}*iAk];
        Kkv = [Kkv;Kk(:,s)];
        Gcheck = Gcheck + Gk{s,k}; % double-check Lemma 2, sum should be an identity matrix
    end
    % update
    Qkb = Qkv * Qk * Qkv' + Kkv * Kkv' .* kron(Rk,ones(n,n)); % matrix Q bold in Lemma 3
    Pikb = kron(eye(m),Pik);
    Wkb{1,k} = Pikb * Wkb{1,k-1} * Pikb' + Qkb; % matrix W bold in Lemma 4
    iWkb = (Wkb{1,k})\eye(n*m); % inverse of matrix W bold
    assign(xlsv,Ak*xls(:,k-1)); % initialize optimizer 
    assign(mu,zeros(n*m,1)); % initialize optimizer 
    assign(var,zeros(n*m,1)); % initialize optimizer 
    varw = diag(kron(phi(:,k),ones(n,1))) * var; % initialize optimizer 
    obj = 0.5 * mu' * iWkb * mu + gammaLS * norm(varw,1); % objective function for the secure least square in Section IV.B, eq. (20)
    F = []; % initialize the constraints
    F = [F, zetav == Gkv * xlsv + mu + varw]; % constraint of the secure least square in Section IV.B, eq. (20)
    opt1 = optimize(F,obj,options); % solve the secure leasqure
    xls(:,k) = value(xlsv); % take the solution of the least square as the state estimate
    muls(:,k) = value(mu); % take the solution of the least square
    varls(:,k) = value(var); % take the solution of the least square
    rmse_xls(k) = norm((xls(:,k) - x(:,k))); % root mean square error of the state estimate with the true state by solving the secure least square 
    for s = 1:m
        if phi(s,k) == 1 && abs(y(s,k) - C(s,:) * xls(:,k)) > 1
            zeta{s}(:,k) = Gk{s,k} * xls(:,k); 
            disp("Attack detected at:" + num2str(s)); % detect anomaly at some sensor
        end
    end
end


%% figure
bus = [2,3,4,5];
figure(1)
set(gcf,'position',[2000,100,1000,600]);
allplots = [];
for ibus = 1:length(bus)
    subplot(1,4,ibus);
    plot(ts,x(bus(ibus),:),'-','LineWidth',3,'Color','black');
    hold on
    plot(ts,xhc(bus(ibus),:),'--','LineWidth',3,'Color','blue');
    hold on
    plot(ts,xh(bus(ibus),:),':','LineWidth',4,'Color','red');
    hold on
    plot(ts,xls(bus(ibus),:),':','LineWidth',4,'Color','#EDB120');
    grid;
    title("Bus " + num2str(bus(ibus)),'FontSize',14); 
    set(gca,'FontSize',16);
end
legend('True state', 'Est. by KF w/o attacks' ,'Est. by KF under attacks',...
    'Est. by (20) with $\gamma = 2$','FontSize',16,'Orientation', 'horizontal','Interpreter','latex');


figure(2)
set(gcf,'position',[2000,400,800,300]);
plot(ts(2:end),rmse_xhc(2:end),'--','LineWidth',4,'Color','blue');
hold on
plot(ts(2:end),rmse_xh(2:end),':','LineWidth',4,'Color','red');
hold on
plot(ts(2:end),rmse_xls(2:end),':','LineWidth',4,'Color','#EDB120');
xlabel('Time (s)');
ylabel('$$||\hat{x}[k] - x[k]||_2$$','Interpreter','latex');
title('Estimation error comparison','FontSize',16);
set(gca,'FontSize',16);
legend('Est. by KF w/o attacks','Est. by KF under attacks','Est. by (20) with $\gamma = 2$ under attacks','FontSize',16,'Interpreter','latex');
grid;


figure(3)
set(gcf,'position',[2000,400,800,300]);
for i = 1:4
    color = rand(1,3);
    plot(att_activ{i,1},i*ones(1,length(att_activ{i,1})),'s','Markersize',8,'MarkerFaceColor', color, 'MarkerEdgeColor', color);
    hold on
end
xlim([0 max(ts)*1.1]);
xlabel('Time (s)');
ylim([0.1 6.2]);
yticks([1, 2, 3, 4]);
yticklabels({'Bus 3', 'Bus 4', 'Bus 5', 'Bus 2'});
set(gca,'FontSize',16);
title('Launched attack samples');
legend('False-data injection','Time-stamp manipulation','Denial-of-service','Fake-data generation','FontSize',16,'NumColumns', 2);
grid;

