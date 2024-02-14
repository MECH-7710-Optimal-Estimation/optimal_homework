%% OPTIMAL - HW0
clear; close all; clc;

% System Variables
J = 10;                         % kg*m^2
b = 1;                          % Nms/rad
fn_o = 50;                      % Hz
wn_o = fn_o*2*pi;               % rad/s
fn_c = 10;                      % Hz
wn_c = fn_c*2*pi;               % rad/s
zeta = 0.7;
fs = 1e3;                       % Hz
Ts = 1/fs;                      % s

% State Space Matrices
A = [0    1;
     0 -b/J];
B = [  0;
     1/J];
C = [1 0];
D = 0;

% Simulation Variables
dt = 1e-3;                      % s
tend = 0.25;
time_c = 0:dt:tend;             % s
time_d = 0:Ts:tend;             % s
u = ones(length(time_c),1);     % Nm

%% QUESTION 1
syms theta(t) T
dTheta = diff(theta);
ddTheta = diff(dTheta);

fprintf('1a) The differential equation of the system is: \n');
EOM = J*ddTheta == (-b*dTheta + T)/J
fprintf('1b) The statespace model of the system is: \n');
sys_ss = ss(A, B, C, D)
[b,a] = ss2tf(A, B, C, D);
Gsys = tf(b,a);
eigenvalues = eigs(A);
fprintf('\n1c) The eigenvalues of the system are [%0.3g %0.3g]\n\n', eigenvalues);

%% QUESTION 2
O = [  C;
     C*A];                  % Observability Matrix
% Check if the observability matrix is full rank denoting
% if the system is observable or not
is_observable = rank(O) == size(A,1);
fprintf("2a) This system is observable: %s\n", string(is_observable));
% Solve the characteristic equation
eigO = roots([1 2*zeta*wn_o wn_o^2]);
% Craft the L matrix.  Work on paper
L = place(A', C', eigO)';
fprintf("2b) The L matrix for the system is [%0.3g %0.3g]'\n\n", L);

% Simulate the system and the estimator.
% Estimates are multiplied by 2 to show the estimator dynamics.
x = ones(2,length(time_c));
x_hat = 2*ones(2,length(time_c));
for i = 2:length(time_c)
    x_hat(:,i) = x_hat(:,i-1) ...
        + ((A - L*C)*x_hat(:,i-1) + B*u(i-1) + L*C*x(:,i-1))*dt;
    x(:,i) = x(:,i-1) + (A*x(:,i-1) + B*u(i-1))*dt;
end

figure();
tcl = tiledlayout(2,1);
title(tcl, '2c) System vs. Estimator');
nexttile();
hold("on");
title("Theta");
plot(time_c,x(1,:));
plot(time_c,x_hat(1,:), '--r');
xlabel("Time (s)");
ylabel("Theta (rad)");
legend('System', 'Estimation');
nexttile();
title("Theta Dot");
hold("on");
plot(time_c,x(2,:));
plot(time_c,x_hat(2,:), '--r');
xlabel("Time (s)");
ylabel("Theta Dot (rad/s)");
legend('System', 'Estimation');

%% QUESTION 3
Cont = [B A*B];                 % Controllability Matrix
% Check if the controllability matrix is full rank denoting
% if the system is controllable or not
is_controllable = rank(Cont) == size(A,1);
fprintf("3a) This system is controllable: %s\n", string(is_controllable));
% Solve the characteristic equation
eigC = roots([1 2*zeta*wn_c wn_c^2]);
% Craft the K matrix.  Work on paper
K = place(A, B, eigC);
fprintf("3b) The K matrix for the system is [%0.3g %0.3g]\n\n", K);

% Simulate the system and the combined controller/estimator.
% Estimates are multiplied by 2 to show the estimator dynamics
x = ones(2,length(time_c));
x_hat = 2*ones(2,length(time_c));
for i = 2:length(time_c)
    x_hat(:,i) = x_hat(:,i-1) ...
        + ((A - B*K - L*C)*x_hat(:,i-1) - B*K*x_hat(:,i-1) + L*C*x(:,i-1))*dt;
    x(:,i) = x(:,i-1) + (A*x(:,i-1) - B*K*x_hat(:,i-1))*dt;
end

figure();
tcl = tiledlayout(2,1);
title(tcl, '3c) System vs. Combined Controller & Estimator');
nexttile();
hold("on");
title("Theta");
plot(time_c, x(1,:));
plot(time_c, x_hat(1,:), '--r');
xlabel("Time (s)");
ylabel("Theta (rad)");
legend('System', 'Estimation');
nexttile();
hold("on");
title("Theta Dot");
plot(time_c, x(2,:));
plot(time_c, x_hat(2,:), '--r');
xlabel("Time (s)");
ylabel("Theta Dot (rad/s)");
legend('System', 'Estimation');

%% QUESTION 4
% Equivalent Compensator in State Space
Acomp = A - B*K - L*C;
Bcomp = L;
Ccomp = -K;
Dcomp = 0;

[b,a] = ss2tf(Acomp, Bcomp, Ccomp, Dcomp);
Gcomp = tf(b,a)
fprintf("4a) The above compensator resembles a lead compensator as it " + ...
    "adds a pole and a zero, as well as damping, to the system\n");
ol_comp = -1*Gsys*Gcomp;                    % Open Loop Compensator & Plant
fprintf("4b) The closed loop transfer function of the system w/ a compensator: \n");
cl_comp = minreal(ol_comp/(1 + ol_comp))    % Closed Loop Compensator & Plant

figure();
opt = bodeoptions;
opt.Title.String = '4c) Closed Loop Compensator Bode';
opt.Title.FontSize = 12;
opt.Title.FontWeight = 'bold';
bode(cl_comp, opt);

% 4d) Gain & Phase Margin
figure();
margin(ol_comp);

%% QUESTION 5
% Convert State-Space system model to discrete
[Ad, Bd, Cd, Dd] = c2dm(A, B, C, D, Ts, 'zoh');
eigenvalues_d = eigs(Ad);                   % Eigs of the discrete system
fprintf('\n5a) The eigenvalues of the discrete system are [%0.3g %0.3g]\n', eigenvalues_d);
[b,a] = ss2tf(Ad, Bd, Cd, Dd);
Gsys_d = tf(b,a);
eigO_d = exp(eigO*Ts);                      % Discrete Observer eigs
eigC_d = exp(eigC*Ts);                      % Discrete Controller eigs
Ld = place(Ad', Cd', eigO_d)';              % Discrete Observer Matrix
fprintf("5b) The L matrix for the discrete system is [%0.3g %0.3g]'\n", Ld);
Kd = place(Ad, Bd, eigC_d);                 % Discrete Controller Matrix
fprintf("5c) The K matrix for the discrete system is [%0.3g %0.3g]'\n", Kd);

% 5d) Discrete Observer/Controller Poles and Zeros
obs_ctr_poles = eigs(Ad - Bd*Kd - Ld*Cd);
fprintf(['5d) The poles of the closed loop ' ...
    'controller/observer system are [%0.4g%+0.4gj %0.4g%+0.4gj]\n'], ...
    real(obs_ctr_poles(1)), imag(obs_ctr_poles(1)), ...
    real(obs_ctr_poles(2)), imag(obs_ctr_poles(2)));

% Discrete Equivalent Compensator in State Space
Acomp_d = Ad - Bd*Kd - Ld*Cd;
Bcomp_d = Ld;
Ccomp_d = -Kd;
Dcomp_d = 0;

[b,a] = ss2tf(Acomp_d, Bcomp_d, Ccomp_d, Dcomp_d);
Gcomp_d = tf(b,a);
ol_comp_d = -1*Gcomp_d*Gsys_d;
fprintf("5e) The closed loop transfer function of the discrete system w/ a compensator: \n\n");
cl_comp_d = minreal(ol_comp_d/(1 + ol_comp_d))

%% QUESTION 6
% Uncompensated Sim
% Continuous
x_c = ones(2,length(time_c)-1);
for i = 2:length(time_c)
    x_c(:,i) = x_c(:,i-1) + (A*x_c(:,i-1) + B*u(i-1))*dt;
end

% Discrete
x_d = ones(2,length(time_d)-1);
for i = 2:length(time_d)
    x_d(:,i) = Ad*x_d(:,i-1) + Bd*u(i-1);
end

figure();
tcl = tiledlayout(2,1);
title(tcl,'Uncompensated Comparison');
nexttile();
hold("on");
title("Continuous vs. Discrete Theta");
plot(time_c, x_c(1,:));
plot(time_d, x_d(1,:), '--r');
xlabel("Time (s)");
ylabel("Theta (rad)");
legend('Continuous', 'Discrete');
nexttile();
hold("on");
title("Continuous vs. Discrete Theta Dot");
plot(time_c, x_c(2,:));
plot(time_d, x_d(2,:), '--r');
xlabel("Time (s)");
ylabel("Theta Dot (rad/s)");
legend('Continuous', 'Discrete');

% Compensator Sim
% Continuous
x_c = ones(2,length(time_c)-1);
x_comp_c = zeros(2,length(time_c)-1);
for i = 2:length(time_c)
    x_c(:,i) = x_c(:,i-1) + (A*x_c(:,i-1) ...
        + B*Ccomp*x_comp_c(:,i-1))*dt;
    x_comp_c(:,i) = x_comp_c(:,i-1) ...
        + (Acomp*x_comp_c(:,i-1) + Bcomp*C*x_c(:,i-1))*dt;
end

% Discrete
x_d = ones(2,length(time_d)-1);
x_comp_d = zeros(2,length(time_d)-1);
for i = 2:length(time_d)
    x_d(:,i) = Ad*x_d(:,i-1) + Bd*Ccomp_d*x_comp_d(:,i-1);
    x_comp_d(:,i) = Acomp_d*x_comp_d(:,i-1) + Bcomp_d*Cd*x_d(:,i-1);
end

figure();
tcl = tiledlayout(2,1);
title(tcl,'Compensator Comparison');
nexttile();
hold("on");
title("Continuous vs. Discrete Theta");
plot(time_c, x_c(1,:));
plot(time_d, x_d(1,:), '--r');
xlabel("Time (s)");
ylabel("Theta (rad)");
legend('Continuous', 'Discrete');
nexttile();
hold("on");
title("Continuous vs. Discrete Theta Dot");
plot(time_c, x_c(2,:));
plot(time_d, x_d(2,:), '--r');
xlabel("Time (s)");
ylabel("Theta Dot (rad/s)");
legend('Continuous', 'Discrete');

fprintf(['\nThe discrete compensator has a small amount of lag as compared\n' ...
    'to the continuous compensator.  The uncompensated systems look to be\n' ...
    'identical.\n'])