%% PART IV
clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/"));

dt = 0.01;
time = 0:dt:10;

A = [-2.62 12;
     -0.96 -2];
Ad = (eye(size(A)) + A.*dt);
B = [14;
      1];
C = [1 0];

sigmaV = sqrt(0.1);

R = sigmaV;
Q = diag([1 0.03]);

del = heaviside(time);

% (A)
% Simulation
x_true1 = zeros(2,length(time));
y1 = zeros(1,length(time));
for i = 2:length(time)
    x_true1(:,i) = x_true1(:, i-1) + (A*x_true1(:, i-1) + B*del(i-1))*dt;
    y1(:,i) = C*x_true1(:,i) + sigmaV*randn;
end

figure();
hold('on');
plot(time, y1, 'LineWidth', 2);
plot(time, x_true1(1,:), 'LineWidth', 2);
title('Simulated Measurements vs. Time');
xlabel('Time (s)');
ylabel('Yaw-Rate Measurement (rad/s)');
legend('Measurement', 'Truth');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p4a_meas.png", 'Resolution', 300);

% Kalman Filter
x1 = zeros(2,length(time));
P = eye(2);
for i = 2:length(time)
    % Time Update
    xp = Ad*x1(:,i-1);
    Pp = Ad*P*Ad' + Q;
    % Kalman Gain
    L = (Pp*C')/(C*Pp*C' + R);
    % Measurement Update
    x1(:,i) = xp + L*(y1(i) - C*xp);
    P = (eye(2) - L*C)*Pp;
end

mean_error = mean(x_true1 - x1,2);
fprintf('Mean Error of the Yaw Rate Estimate: %0.5g\n', mean_error(1));
fprintf('Mean Error of the Slip Angle Estimate: %0.5g\n\n', mean_error(2));

figure();
t = tiledlayout(2,1);
nexttile();
hold('on');
plot(time, x1(1,:), '--', 'LineWidth', 2);
plot(time, x_true1(1,:), 'LineWidth', 2);
title('Yaw Rate vs. Time');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

nexttile();
hold('on');
plot(time, x1(2,:), '--', 'LineWidth', 2);
plot(time, x_true1(2,:), 'LineWidth', 2);
title('Side Slip Angle vs. Time');
xlabel('Time (s)');
ylabel('Side Slip Angle (rad)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p4a_kf.png", 'Resolution', 300);

%% (B)
A2 = [-2.42  4;
      -0.99 -2];
Ad2 = (eye(2) + A*dt);
B2 = [18
       1];
C2 = C;

Q2 = diag([2 0.06]);

% Simulation
x_true2 = zeros(2,length(time));
y2 = zeros(1,length(time));
for i = 2:length(time)
    x_true2(:,i) = x_true2(:, i-1) + (A*x_true2(:, i-1) + B*del(i-1)).*dt;
    y2(:,i) = C*x_true2(:,i) + sigmaV*randn;
end

% Kalman Filter
x2 = zeros(2,length(time));
P = eye(2);
for i = 2:length(time)
    % Time Update
    xp = Ad*x2(:,i-1);
    Pp = Ad*P*Ad' + Q2;
    % Kalman Gain
    L = (Pp*C')/(C*Pp*C' + R);
    % Measurement Update
    x2(:,i) = xp + L*(y2(i) - C*xp);
    P = (eye(2) - L*C)*Pp;
end

mean_error = mean(x_true2 - x2, 2);
fprintf('Mean Error of the Yaw Rate Estimate: %0.5g\n', mean_error(1));
fprintf('Mean Error of the Slip Angle Estimate: %0.5g\n\n', mean_error(2));

figure();
tiledlayout(2,1);
nexttile();
hold('on');
plot(time, x2(1,:), '--', 'LineWidth', 2);
plot(time, x_true2(1,:), 'LineWidth', 2);
title('Yaw Rate vs. Time');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

nexttile();
hold('on');
plot(time, x2(2,:), '--', 'LineWidth', 2);
plot(time, x_true2(2,:), 'LineWidth', 2);
title('Side Slip Angle vs. Time');
xlabel('Time (s)');
ylabel('Side Slip Angle (rad)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p4b_kf.png", 'Resolution', 300);

%% (C)
sigmaN = sqrt(0.5);
R3 = diag([sigmaV sigmaN]);

% (D)
A3 = A;
B3 = B;
C3 = eye(2);

Q3 = diag([2, 0.02]);

% Simulation
x_true3 = zeros(2,length(time));
y3 = zeros(2,length(time));
for i = 2:length(time)
    x_true3(:,i) = x_true3(:, i-1) + (A*x_true3(:, i-1) + B*del(i-1)).*dt;
    y3(:,i) = C3*x_true3(:,i) + [sigmaV sigmaN]*randn(2,1);
end

% Kalman Filter
x3 = zeros(2,length(time));
P = eye(2);
for i = 2:length(time)
    % Time Update
    xp = Ad*x3(:,i-1);
    Pp = Ad*P*Ad' + Q3;
    % Kalman Gain
    L = (Pp*C3')/(C3*Pp*C3' + R3);
    % Measurement Update
    x3(:,i) = xp + L*(y3(:,i) - C3*xp);
    P = (eye(2) - L*C3)*Pp;
end

mean_error = mean(x_true3 - x3, 2);
fprintf('Mean Error of the Yaw Rate Estimate: %0.5g\n', mean_error(1));
fprintf('Mean Error of the Slip Angle Estimate: %0.5g\n\n', mean_error(2));

snrratio_meas = snr(x_true3(2,:), y3(2,:));
snratio_estimate = snr(y3(2,:), x3(2,:));

figure();
tiledlayout(2,1);
nexttile();
hold('on');
plot(time, x3(1,:), '--', 'LineWidth', 2);
plot(time, x_true3(1,:), 'LineWidth', 2);
title('Yaw Rate vs. Time');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

nexttile();
hold('on');
plot(time, x3(2,:), '--', 'LineWidth', 2);
plot(time, x_true3(2,:), 'LineWidth', 2);
title('Side Slip Angle vs. Time');
xlabel('Time (s)');
ylabel('Side Slip Angle (rad)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p4d_kf.png", 'Resolution', 300);