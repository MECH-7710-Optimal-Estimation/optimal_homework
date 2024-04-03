clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/"));

data = readtable("data/hw3_2.txt");
t = data.Var1;
dt = mean(diff(t));
y1 = data.Var2;

R = 1;
Qd1 = 0;

Ac = 0;
Ad = eye(size(Ac)) - Ac*dt;
C = 1;
Bw = 0;

x1 = zeros(1,length(y1));
P1 = 1;
L1 = zeros(1,length(t));
for i = 2:length(t)
    % Time Update
    xp = Ad*x1(:,i-1);
    Pp = Ad*P1*Ad' + Qd1;
    % Kalman Gain
    L1(:,i) = (Pp*C)/(C'*Pp*C + R);
    % Measurement Update
    x1(:,i) = xp + L1(:,i)*(y1(i) - C*xp);
    P1 = (eye(size(Ad)) - L1(:,i)*C)*Pp;
end

figure();
hold('on');
title('Bias & Kalman Gain vs. Time');
plot(t, x1, 'LineWidth', 2);
plot(t, L1, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bias & Kalman Gain');
legend('Bias', 'Kalman Gain')
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2a.png", 'Resolution', 300);

Qd2 = 0.0001;
Qd3 = 0.01;

P2 = 1;
x2 = zeros(1,length(y1));
L2 = zeros(1,length(t));
P3 = 1;
x3 = zeros(1,length(y1));
L3 = zeros(1,length(t));
for i = 2:length(t)
    % Time Update
    xp2 = Ad*x2(:,i-1);
    Pp2 = Ad*P2*Ad' + Qd2;
    xp3 = Ad*x3(:,i-1);
    Pp3 = Ad*P3*Ad' + Qd3;
    % Kalman Gain
    L2(:,i) = (Pp2*C)/(C'*Pp2*C + R);
    L3(:,i) = (Pp3*C)/(C'*Pp3*C + R);
    % Measurement Update
    x2(:,i) = xp2 + L2(:,i)*(y1(i) - C*xp2);
    P2 = (eye(size(Ad)) - L2(:,i)*C)*Pp2;
    x3(:,i) = xp3 + L3(:,i)*(y1(i) - C*xp3);
    P3 = (eye(size(Ad)) - L3(:,i)*C)*Pp3;
end

figure();
hold('on');
title('Bias & Kalman Gain vs. Time');
plot(t, x2, 'LineWidth', 2);
plot(t, L2, 'LineWidth', 2);
plot(t, x3, 'LineWidth', 2);
plot(t, L3, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bias & Kalman Gain');
legend('Bias; Q_D = 0.0001', 'L_K; Q_D = 0.0001', ...
       'Bias; Q_D = 0.01', 'L_K; Q_D = 0.01')
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2b.png", 'Resolution', 300);

% Tradeoff: Settle Time to Noise

numd2 = sqrt(Qd2);
dend2 = [1 (-1+sqrt(Qd2))];
y_filt2 = filter(numd2, dend2, y1, y1(1));
numd3 = sqrt(Qd3);
dend3 = [1 (-1+sqrt(Qd3))];
y_filt3 = filter(numd3, dend3, y1, y1(1));

figure();
hold('on');
title('Bias & Filter vs. Time');
plot(t, x2, 'LineWidth', 2);
plot(t, y_filt2, '--', 'LineWidth', 2);
plot(t, x3, 'LineWidth', 2);
plot(t, y_filt3, '--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bias & Filter');
legend('Bias; Q_D = 0.0001', 'y_{filt}; Q_D = 0.0001', ...
       'Bias; Q_D = 0.01', 'y_{filt}; Q_D = 0.01');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2c.png", 'Resolution', 300);