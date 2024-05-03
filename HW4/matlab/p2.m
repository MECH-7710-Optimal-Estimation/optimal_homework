clear; close all; clc;

N = 1000;                       % Number of Samples

numd = 0.25*[1 -0.8];           % Discrete TF Numerator
dend = [1 -1.9 0.95];           % Discrete TF Denominator
tfd = tf(numd, dend, -1);       % Discrete TF
u = randn(N,1);                 % Input Noise
y = dlsim(numd,dend,u);         % Output Simulation
sigma1 = 0.01;                  % Output Noise Standard Deviation
Y = y + sigma1*randn(N,1);      % Output w/ Noise

H = [-Y(2:end-1) -Y(1:end-2) u(2:end-1) u(1:end-2)];
x = pinv(H)*Y(3:end);
numLS = x(3:4)';
denLS = [1 x(1:2)'];
tfLS = tf(numLS, denLS, -1);
yLS = dlsim(numLS, denLS, u);

tfARX = arx(u, Y, [2,3,0]);

figure();
bode(tfLS, '-r');
hold('on');
bode(tfARX, '--');
bode(tfd, '-.k');
title('Truth v. SysID Bodes');
legend('LS', 'ARX', 'Simulated');

sigma1a = 1;
Ya = y + sigma1a*randn(N,1);      % Output w/ Noise
tfARXa = arx(u, Ya, [2,3,0]);

figure();
bode(tfARXa, '-r');
hold('on');
bode(tfd, '-.k');
title('Truth v. ARX Bode; 1\sigma=1 Noise');
legend('ARX', 'Simulated');

tfARXa1 = arx(u, Ya, [3,4,0]);

figure();
bode(tfARXa1, '-r');
hold('on');
bode(tfd, '-.k');
title('Truth v. SysID Bodes; A: 3rd Order; B: 4th Order');
legend('ARX', 'Simulated');

tfARXa2 = arx(u, Ya, [10,10,0]);

figure();
bode(tfARXa2, '-r');
hold('on');
bode(tfd, '-.k');
title('Truth v. SysID Bodes; A: 3rd Order; B: 4th Order');
legend('ARX', 'Simulated');

tfARXa3 = arx(u, Ya, [20,20,0]);

figure();
bode(tfARXa3, '-r');
hold('on');
bode(tfd, '-.k');
title('Truth v. SysID Bodes; A: 3rd Order; B: 4th Order');
legend('ARX', 'Simulated');

% Could not identify model after high noise.

%% (B)
tfARMAX = armax(u,Ya,[2,3,2,0]);
tfBJ = bj(u,y,[3,3,3,2,0]);
tfIV4 = iv4(u, y, [2,3,0]);

figure();
bode(tfARMAX, '-or');
hold('on');
bode(tfBJ, '-*b');
bode(tfd, '-squarek');
bode(tfd, '-xk');
title('Truth v. SysID Bodes');
legend('ARMAX', 'Box-Jenkins', 'IV4', 'Simulated');

% All three function better than ARX