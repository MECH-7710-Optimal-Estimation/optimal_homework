%% PART 3
clear; close all; clc;

N1 = 10;                    % Number of Samples
Nmc = 1000;                 % Number of Monte Carlos
sigma_n = 0.3;              % Standard Deviation [deg/s]
t1 = 0:(N1-1);              % Time [s]
w1 = 2*pi/N1;               % Frequency [rad/s]
a = 1;                      % Gyroscope Bias Factor
b = 0;                      % Gyroscope Bias Factor

% Simulation
n1 = sigma_n*randn(N1,1);       % Noise
r1 = 100*sin(w1*t1)';           % Sine Wave
g1 = a*r1 + b*ones(N1,1) + n1;  % Gyroscope Measurements

% (A)
% Least Squares
Ha = [r1 ones(N1,1)];        % Geometry Matrix
xa = pinv(Ha)*g1;            % State Estimate
fprintf('a) 10 Sample; Individual Least Squares Estimate: [%0.3g %0.3g]\n', xa);

% (B)
xb = zeros(2,Nmc);
for i = 1:Nmc
    n1 = sigma_n*randn(N1,1);       % Noise
    g1 = a*r1 + b*ones(N1,1) + n1;  % Gyroscope Measurements
    Hb = [r1 ones(N1,1)];           % Geometry Matrix
    xb(:,i) = pinv(Hb)*g1;          % State Estimate
end
sigmab = std(xb,0,2);       % Sample Sigma
meanb = mean(xb,2);         % Sample Mean
sigmab_an = sigma_n/sqrt(N1);
fprintf('b) 10 Sample; Monte Carlo Least Squares Mean: [%0.3g %0.3g]\n', meanb);
fprintf('b) 10 Sample; Monte Carlo Least Squares Sigma: [%0.3g %0.3g]\n', sigmab);

N2 = 1000;
t2 = 0:(N2-1);              % Time [s]
w2 = 2*pi/N2;               % Frequency [Hz]

% Simulation
n2 = sigma_n*randn(N2,1);       % Noise
r2 = 100*sin(w2*t2)';           % Sine Wave
g2 = a*r2 + b*ones(N2,1) + n2;  % Gyroscope Measurements

% (C)
% Least Squares
Hc = [r2 ones(N2,1)];        % Geometry Matrix
xc = pinv(Hc)*g2;            % State Estimate
fprintf('c) 1000 Sample; Individual Least Squares Estimate: [%0.3g %0.3g]\n', xc);

xc = zeros(2,Nmc);
for i = 1:Nmc
    n2 = sigma_n*randn(N2,1);       % Noise
    g2 = a*r2 + b*ones(N2,1) + n2;  % Gyroscope Measurements
    Hc = [r2 ones(N2,1)];           % Geometry Matrix
    xc(:,i) = pinv(Hc)*g2;          % State Estimate
end
sigmac = std(xc,0,2);       % Sample Sigma
meanc = mean(xc,2);         % Sample Mean
fprintf('c) 1000 Sample; Monte Carlo Least Squares Mean: [%0.3g %0.3g]\n', meanc);
fprintf('c) 1000 Sample; Monte Carlo Least Squares Sigma: [%0.3g %0.3g]\n', sigmac);

% (D)
xr = zeros(2,N2);

% Simulation
nd = sigma_n*randn(N2,1);       % Noise
rd = 100*sin(w2*t2)';           % Sine Wave
gd = a*rd + b*ones(N2,1) + nd;  % Gyroscope Measurements

R = diag([sigma_n, sigma_n]);
batch = 2;

Hd = Hc(1:batch,:);
xr(:,1) = pinv(Hd)*gd(1:batch);
P = inv(Hd'*Hd);

idx = 2;
for i = batch:batch:N2-batch
    Hd = Hc(i:i+batch-1,:);
    % Kk = P*Hk'*inv(Hk*P*Hk' + R);
    % P = (eye(2) - Kk*Hk)*P;
    P = inv(inv(P) + Hd'*inv(R)*Hd);
    K = P*Hd'*inv(R);
    xr(:,idx) = xr(:,idx-1) + K*(gd(i:i+batch-1) - Hd*xr(:,idx-1));
end

figure();
plot(xr');

%% PART 4
clear; close all; clc;

N = 1000;                       % Number of Samples

numd = 0.25*[1 -0.8];           % Discrete TF Numerator
dend = [1 -1.9 0.95];           % Discrete TF Denominator
tfd = tf(numd, dend);           % Discrete TF
u = randn(N,1);                 % Input Noise
y = dlsim(numd,dend,u);         % Output Simulation
sigma1 = 0.01;                  % Output Noise Standard Deviation
Y = y + sigma1*randn(N,1);      % Output w/ Noise

Ha = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
xa = pinv(Ha)*Y(3:end);

numID = xa(3:4)';
denID = [1 xa(1:2)'];
tfID = tf(numID, denID);
yID = dlsim(numID,denID,u);

figure();
bode(tfd, '-b');
hold('on');
bode(tfID, '--r')
title('Truth v. SysID Bode');
legend('Simulated', 'SysID');

figure();
hold('on');
plot(y, 'LineWidth', 2);
plot(yID, '--', 'LineWidth', 2);
title('Simulated Truth v. Simulated SysID');
legend('Truth', 'SysID');

N_mc = 10;
x_mc1 = zeros(4,N_mc);
for i = 1:N_mc
    u = randn(N,1);
    y = dlsim(numd,dend,u);
    Y = y + sigma1*randn(N,1);
    Ha = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
    x_mc1(:,i) = pinv(Ha)*Y(3:end);
end
sigma_x_an1 = diag(sigma1*inv(Ha'*Ha));

sigma_x_mc1 = std(x_mc1,1,2);
mean_x_mc1 = mean(x_mc1,2);

sigma1 = 0.9;
x_mc2 = zeros(4,N_mc);
for i = 1:N_mc
    u = randn(N,1);
    y = dlsim(numd,dend,u);
    Y = y + sigma1*randn(N,1);
    Ha = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
    x_mc2(:,i) = pinv(Ha)*Y(3:end);
end
sigma_x_an2 = diag(sigma1*inv(Ha'*Ha));

sigma_x_mc2 = std(x_mc2,1,2);
mean_x_mc2 = mean(x_mc2,2);

%% PART V
clear;
Tw = 1;
A = 0.1;
wc = 2;     % [rad/s]
w1 = -10:0.01:10; % [rad/s]

syms w
Sw = piecewise(abs(w) <= wc, A, abs(w) > wc, 0);
G = abs(1./(Tw*w1 + 1));

for i = 1:length(w1)
    if abs(w1(i)) < wc
        Sy1(i) = (A)/(Tw*w1(i) + 1);
    else
        Sy1(i) = 0;
    end
end

Sy2 = A./(Tw.*w1 + 1);

figure();
fplot(Sw, [-10 10], 'LineWidth', 2);
xlabel('Frequency [rad/s]');
ylabel('Magnitude');
title('PSD of Noise');
xline(-wc, '--k','LineWidth', 2);
xline(wc, '--k','LineWidth', 2);
ylim([-A*2 A*2])
ax = gca;
ax.FontSize = 18;

figure();
plot(w1, G, 'LineWidth', 2);
xlabel('Frequency [rad/s]');
ylabel('|G(jw)|');
title('|G(jw)| vs. Frequency');
ax = gca;
ax.FontSize = 18;

figure();
hold('on');
plot(w1,abs(Sy1), 'LineWidth', 2);
plot(w1,abs(Sy2), '--', 'LineWidth', 2);
xline(-wc, '--k','LineWidth', 2);
xline(wc, '--k','LineWidth', 2);
xlabel('Frequency [rad/s]');
ylabel('Sy');
title('Sy vs. Frequency');
legend('Sy1', 'Sy2');
ax = gca;
ax.FontSize = 18;