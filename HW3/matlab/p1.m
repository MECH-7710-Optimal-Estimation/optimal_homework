clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/"));

Ts = 0.1;           % Sample Period [s]
fs = 1/Ts;          % Sample Rate [Hz]

Acl = [ 0    1;
       -1 -1.4];    % Continuous State Transtion Matrix
C = [1 0];          % Observation Matrix
Bw = [0 1];         % Input Noise Matrix

sigmaV = 1;         % Sensor Noise Sigma
sigmaQ = 2;         % Process Noise Sigma

% (A)
time = 0:Ts:100;
x_true1 = zeros(2,length(time));
y1 = zeros(2,length(time));
for i = 2:length(time)
    x_true1(:,i) = x_true1(:,i-1) + (Acl*x_true1(:,i-1) + Bw*sigmaQ*randn(2,1))*Ts;
    y1(:,i) = C*x_true1(:,i) + sigmaV*randn(2,1);
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold('on');
title('A) Simulated Measurement vs. Time');
plot(time, y1, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Measurement');
legend('Position [m]', 'Velocity [m/s]')
ax = gca;
ax.FontSize = 18;

% exportgraphics(gcf, currentFolder + "/../figures/p1a.png", 'Resolution', 300);

% (B)
R = sigmaV^2.*eye(size(Acl));   % Measurement Covariance Matrix
Q = sigmaQ^2.*eye(size(Acl));   % Continuous Process Covariance Matrix
Qd1 = Bw*Q*Bw'*Ts;              % Discrete Process Covariance Matrix
Ad = eye(size(Acl)) + Acl.*Ts;  % Discrete State Transition Matrix

%% (C)
x1 = 10.*ones(2,length(time));  % State Estimates
L1 = zeros(2,length(time));     % Kalman Gain
P1 = zeros(2,2,length(time));   % Covariance Matrix
P1(:,:,1) = eye(2);
Pp1 = zeros(2,2,length(time));   % A Priori Covariance Matrix
for i = 2:length(time)
    % Time Update
    xp = Ad*x1(:,i-1);
    Pp1(:,:,i) = Ad*P1(:,:,i-1)*Ad' + Qd1;
    % Kalman Gain
    L1(:,i) = (Pp1(:,:,i)*C')/(C*Pp1(:,:,i)*C' + sigmaV^2);
    % Measurement Update
    x1(:,i) = xp + L1(:,i)'*(y1(:,i-1) - C*xp);
    P1(:,:,i) = (eye(size(Pp1(:,:,i))) - L1(:,i)*C)*Pp1(:,:,i);
end

figure();
tiledlayout(2,1);
nexttile();
hold('on');
plot(time, x1(1,:), '--', 'LineWidth', 2);
plot(time, x_true1(1,:), 'LineWidth', 2);
title('Position vs. Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

nexttile();
hold('on');
plot(time, x1(2,:), '--', 'LineWidth', 2);
plot(time, x_true1(2,:), 'LineWidth', 2);
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

% exportgraphics(gcf, currentFolder + "/../figures/p1c1.png", 'Resolution', 300);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold('on');
title('C) Kalman Gain vs. Time');
plot(time, L1, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Kalman Gain');
legend('Gain on Position', 'Gain on Velocity');
ax = gca;
ax.FontSize = 18;

% exportgraphics(gcf, currentFolder + "/../figures/p1c2.png", 'Resolution', 300);

%% (D)
Ppss = squeeze(Pp1(:,:,end));
Pss = squeeze(P1(:,:,end));     % Steady State Covariance Matrix
Lss = L1(:,end);                % Steady State Kalman Gain
x2 = 10.*ones(2,length(time));  % State Estimate
for i = 2:length(time)
    % Time Update
    xp = Ad*x2(:,i-1);
    % Measurement Update
    x2(:,i) = xp + Lss'*(y1(:,i-1) - C*xp);
end

poles = eigs(Ad - Lss*C);

figure();
tiledlayout(2,1);
nexttile();
hold('on');
plot(time, x2(1,:), '--', 'LineWidth', 2);
plot(time, x_true1(1,:), 'LineWidth', 2);
title('Position vs. Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

nexttile();
hold('on');
plot(time, x2(2,:), '--', 'LineWidth', 2);
plot(time, x_true1(2,:), 'LineWidth', 2);
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Estimate', 'Truth');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p1d.png", 'Resolution', 300);

N1 = norm(std(x_true1 - x1, 0, 2));     % Normal Kalman Filter
N2 = norm(std(x_true1 - x2, 0, 2));     % SS Kalman Filter

%% (E)
dQR = 0.5;
Q_range = dQR:dQR:5;
R_range = dQR:dQR:5;
Qidx = 1;
Ridx = 1;
mean_error = zeros(length(Q_range), length(R_range));
for Q = Q_range
    for R = R_range
        Lss = dlqr(Ad, Bw', Q, R);
        x3 = zeros(2,length(time));
        for i = 2:length(time)
            % Time Update
            xp = Ad*x3(:,i-1);
            % Measurement Update
            x3(:,i) = xp + Lss*(y1(:,i-1) - C*xp);
        end
        % mean_error(Qidx, Ridx) = norm(mean(x_true1 - x3, 2));
        % mean_error(Qidx, Ridx) = norm(std(x_true1 - x3, 0, 2));
        mean_error(Qidx, Ridx) = sqrt(std(x_true1(1,:).^2 - x3(1,:).^2) + std(x_true1(2,:).^2 - x3(2,:).^2));
        Ridx = Ridx + 1;
    end
    Qidx = Qidx + 1;
    Ridx = 1;
end

minError = min(mean_error, [], "all");
[minQidx, minRidx] = find(mean_error == minError);

minQ = Q_range(minQidx);
minR = R_range(minRidx);