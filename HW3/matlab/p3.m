clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/"));
data = readtable("data/hw3_3.txt");

time = data.Var1;
meas_east = data.Var2;
meas_north = data.Var3;
meas_psi = data.Var4;
meas_psi_dot = data.Var5;
meas_vel = data.Var6;

y = [meas_east, meas_north, meas_psi, meas_psi_dot, meas_vel];
dt = mean(diff(time));
N = length(time);
fs = 1/dt;

% (A)
xKF = zeros(5,length(time));
P = eye(5);
Qd = diag([0.01 0.01 0.01 0.001 0.001]);
Rd = diag([0.1 0.1 0.1 0.1 0.1]);

for i = 2:length(time)
    A = zeros(5);
    A(1,3) =  (y(i-1,5) - xKF(4,i-1))*cos(xKF(3,i-1));
    A(1,4) = -sin(xKF(3,i-1));
    A(2,3) = -(y(i-1,5) - xKF(4,i-1))*sin(xKF(3,i-1));
    A(2,4) = -cos(xKF(3,i-1));
    A(3,5) = -1;
    Phi = eye(5) - A*dt;
    % Time Update
    xp = Phi*xKF(:,i-1);
    Pp = Phi*P*Phi' + Qd;
    % Kalman Gain
    H = eye(5);
    H(4,4) = 0;
    H(4,5) = 1;
    H(5,4) = 1;
    H(5,5) = 0;
    L = (Pp*H')/(H*Pp*H' + Rd);
    % Measurement Update
    xKF(:,i) = xp + L*(y(i,:)' - H*xp);
    P = (eye(5) - L*H)*Pp;
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
tiledlayout(2,2);
nexttile();
plot(time, xKF(1:2,:), 'LineWidth', 2);
title('East & North vs. Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('E', 'N');
ax = gca;
ax.FontSize = 18;

nexttile();
plot(time, rad2deg(xKF(3,:)), 'LineWidth', 2);
title('Heading vs. Time');
xlabel('Time (s)');
ylabel('Heading (deg)');
ax = gca;
ax.FontSize = 18;

nexttile([1,2]);
plot(time, xKF(4:5,:), 'LineWidth', 2);
title('Biases vs. Time');
xlabel('Time (s)');
ylabel('Bias (m/s | rad/s)');
legend('b_r', 'b_g');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p3a.png", 'Resolution', 300);

%% (B)
batch = 1;
xRLS = zeros(5,N);
H = eye(5);
H(4,4) = 0;
H(4,5) = 1;
H(5,4) = 1;
H(5,5) = 0;

xRLS(:,1) = pinv(H)*y(1:batch,:)';
Qd2 = diag([0.01 0.01 0.01 0.001 0.001]);
Rd2 = diag([0.1 0.1 0.1 0.1 0.1]);
Qd2(4,4) = 0;
Qd2(5,5) = 0;
% P = inv(H'*Qd*H);
P = eye(5);
for i = 2:round((N/batch))
    start = i*batch; stop = start + batch - 1;
    K = (P*H')/(H*P*H' + Rd2);
    P = (eye(5) - K*H)*P;
    xRLS(:,i) = xRLS(:,i-1) + K*(y(start:stop,:)' - H*xRLS(:,i-1));
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
tiledlayout(2,2);
nexttile();
plot(time, xRLS(1:2,:), 'LineWidth', 2);
title('East & North vs. Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('E', 'N');
ax = gca;
ax.FontSize = 18;

nexttile();
plot(time, rad2deg(xRLS(3,:)), 'LineWidth', 2);
title('Heading vs. Time');
xlabel('Time (s)');
ylabel('Heading (deg)');
ax = gca;
ax.FontSize = 18;

nexttile([1,2]);
plot(time, xRLS(4:5,:), 'LineWidth', 2);
title('Biases vs. Time');
xlabel('Time (s)');
ylabel('Bias (m/s | rad/s)');
legend('b_r', 'b_g');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p3b.png", 'Resolution', 300);

%% (C)
xKF2 = zeros(5,length(time));
P = eye(5);
Qd = diag([0.01 0.01 0.01 0.001 0.001]);
Rd = diag([0.1 0.1 0.1 0.1 0.1]);
for i = 2:length(time)
    A = zeros(5);
    A(1,3) =  (y(i-1,5) - xKF2(4,i-1))*cos(xKF2(3,i-1));
    A(1,4) = -sin(xKF2(3,i-1));
    A(2,3) = -(y(i-1,5) - xKF2(4,i-1))*sin(xKF2(3,i-1));
    A(2,4) = -cos(xKF2(3,i-1));
    A(3,5) = -1;
    Phi = eye(5) - A*dt;
    % Time Update
    xp = Phi*xKF2(:,i-1);
    Pp = Phi*P*Phi' + Qd;
    % Kalman Gain
    H = eye(5);
    H(4,4) = 0;
    H(5,5) = 0;
    if i < (N - (40/dt))
        H(4,5) = 1;
        H(5,4) = 1;
    end
    L = (Pp*H')/(H*Pp*H' + Rd);
    % Measurement Update
    xKF2(:,i) = xp + L*(y(i,:)' - H*xp);
    P = (eye(5) - L*H)*Pp;
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
tiledlayout(2,2);
nexttile();
plot(time, xKF2(1:2,:), 'LineWidth', 2);
title('East & North vs. Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('E', 'N');
ax = gca;
ax.FontSize = 18;

nexttile();
plot(time, rad2deg(xKF2(3,:)), 'LineWidth', 2);
title('Heading vs. Time');
xlabel('Time (s)');
ylabel('Heading (deg)');
ax = gca;
ax.FontSize = 18;

nexttile([1,2]);
plot(time, xKF2(4:5,:), 'LineWidth', 2);
title('Biases vs. Time');
xlabel('Time (s)');
ylabel('Bias (m/s | rad/s)');
legend('b_r', 'b_g');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p3c.png", 'Resolution', 300);