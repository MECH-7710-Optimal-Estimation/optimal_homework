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
sigmab_an = sqrt(sigma_n^2*inv(Hb'*Hb));
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

sigmac_an = sqrt(sigma_n^2*inv(Hc'*Hc));
fprintf('c) 1000 Sample; Monte Carlo Least Squares Mean: [%0.3g %0.3g]\n', meanc);
fprintf('c) 1000 Sample; Monte Carlo Least Squares Sigma: [%0.3g %0.3g]\n', sigmac);

% (D)
N3 = 1000;
batch = 2;
t3 = 0:(N3-1);              % Time [s]
w3 = 2*pi/N3;               % Frequency [Hz]
Rk = diag([sigma_n, sigma_n]);

% Simulation
n3 = sigma_n*randn(N3,1);       % Noise
r3 = 100*sin(w3*t3)';           % Sine Wave
g3 = a*r3 + b*ones(N3,1) + n3;  % Gyroscope Measurements

x3 = zeros(2,(N3/batch)-1);
P3 = zeros(2,2,(N3/batch)-1);

H3 = [r3 ones(N3,1)];

x3(:,1) = pinv(H3(1:batch,:))*g3(1:batch);
P3(:,:,1) = inv(H3(1:batch,:)'*H3(1:batch,:));

for i = 2:(N3/batch)-1
    start = i*batch; stop = start + batch - 1;
    Hk = H3(start:stop, :);
    K = P3(:,:,i-1)*Hk'*inv(Hk*P3(:,:,i-1)*Hk' + Rk);
    P3(:,:,i) = (eye(2) - K*Hk)*P3(:,:,i-1);
    x3(:,i) = x3(:,i-1) + K*(g3(start:stop) - Hk*x3(:,i-1));
end

sigma3 = std(x3,0,2);
mean3 = mean(x3,2);

figure();
t = tiledlayout(2,1);
title(t, 'Recursive Least Squares Estimates vs. Time', 'FontSize', 20);
nexttile(t);
hold('on');
plot(x3(1,:), 'LineWidth', 2);
yline(sigma3(1) + mean3(1), '--r', 'LineWidth', 2);
yline(-sigma3(1) + mean3(1), '--r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Scale Factor Estimate');
title('Scale Factor Estimate vs. Time');
ax = gca;
ax.FontSize = 18;

nexttile(t);
plot(x3(2,:), 'LineWidth', 2);
yline(sigma3(2) + mean3(2), '--r', 'LineWidth', 2);
yline(-sigma3(2) + mean3(2), '--r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bias Estimate');
title('Bias Estimate vs. Time');
ax = gca;
ax.FontSize = 18;

leg = legend('State Estimate', '1\sigma', 'Orientation', 'horizontal');
leg.Layout.Tile = 'north';
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, 'figures/p3d_state.png', 'Resolution', 300);

figure();
t = tiledlayout(2,1);
title(t, 'Recursive Least Squares Variances vs. Time', 'FontSize', 20);
nexttile(t);
hold('on');
plot(squeeze(P3(1,1,:)), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Scale Factor Variance');
title('Scale Factor Variance vs. Time');
ax = gca;
ax.FontSize = 18;

nexttile(t);
plot(squeeze(P3(2,2,:)), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bias Variance');
title('Bias Variance vs. Time');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, 'figures/p3d_variance.png', 'Resolution', 300);

%% PART 4
clear; close all; clc;

N = 1000;                       % Number of Samples

numd = 0.25*[1 -0.8];           % Discrete TF Numerator
dend = [1 -1.9 0.95];           % Discrete TF Denominator
tfd = tf(numd, dend, -1);           % Discrete TF
u = randn(N,1);                 % Input Noise
y = dlsim(numd,dend,u);         % Output Simulation
sigma1 = 0.01;                  % Output Noise Standard Deviation
Y = y + sigma1*randn(N,1);      % Output w/ Noise

Ha = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
xa = pinv(Ha)*Y(3:end);

numID = xa(3:4)';
denID = [1 xa(1:2)'];
tfID = tf(numID, denID, -1);
yID = dlsim(numID,denID,u);
snr1 = snr(y, Y-y);

figure();
bode(tfd, '-b');
hold('on');
bode(tfID, '--r')
title('Truth v. SysID Bode');
legend('Simulated', 'SysID');

exportgraphics(gcf, 'figures/p4b_bode.png', 'Resolution', 300);

figure();
hold('on');
plot(y, 'LineWidth', 2);
plot(yID, '--', 'LineWidth', 2);
title('Simulated Truth v. Simulated SysID');
legend('Truth', 'SysID');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, 'figures/p4b_tf.png', 'Resolution', 300);

% snrVal = snr(y, Y - y)

N_mc = 10;
x_mc1 = zeros(4,N_mc);
for i = 1:N_mc
    u = randn(N,1);
    y = dlsim(numd,dend,u);
    Y = y + sigma1*randn(N,1);
    Ha = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
    x_mc1(:,i) = pinv(Ha)*Y(3:end);
end
sigma_x_an1 = sqrt(diag(sigma1.^2*inv(Ha'*Ha)));

sigma_x_mc1 = std(x_mc1,1,2);
mean_x_mc1 = mean(x_mc1,2);

sigma2 = 0.9;
Y = y + sigma2*randn(N,1);      % Output w/ Noise

Hd = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
xd = pinv(Hd)*Y(3:end);

numID = xd(3:4)';
denID = [1 xd(1:2)'];
tfID = tf(numID, denID, -1);
yID = dlsim(numID,denID,u);
snr2 = snr(y, Y-y);

figure();
bode(tfd, '-b');
hold('on');
bode(tfID, '--r')
title('Truth v. SysID Bode');
legend('Simulated', 'SysID');

exportgraphics(gcf, 'figures/p4d_bode.png', 'Resolution', 300);

figure();
hold('on');
plot(y, 'LineWidth', 2);
plot(yID, '--', 'LineWidth', 2);
title('Simulated Truth v. Simulated SysID');
legend('Truth', 'SysID');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, 'figures/p4d_tf.png', 'Resolution', 300);

x_mc2 = zeros(4,N_mc);
for i = 1:N_mc
    u = randn(N,1);
    y = dlsim(numd,dend,u);
    Y = y + sigma2*randn(N,1);
    Hd = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
    x_mc2(:,i) = pinv(Hd)*Y(3:end);
end
sigma_x_an2 = sqrt(diag(sigma2^2*inv(Hd'*Hd)));

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

exportgraphics(gcf, 'figures/p5_noise_psdl.png', 'Resolution', 300);

figure();
plot(w1, G, 'LineWidth', 2);
xlabel('Frequency [rad/s]');
ylabel('|G(jw)|');
title('|G(jw)| vs. Frequency');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, 'figures/p5_gjw.png', 'Resolution', 300);

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

exportgraphics(gcf, 'figures/p5_sy.png', 'Resolution', 300);