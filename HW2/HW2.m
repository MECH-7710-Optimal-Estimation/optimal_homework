%% PART 3
clear; close all; clc;

N = 10;                     % Number of Samples
sigma_n = 0.3;              % Standard Deviation [deg/s]
tf = 1;                     % Final Time [s]
t = linspace(0,tf,N);       % Time [s]
f = 1/N;                    % Frequency [Hz]
w = 2*pi*f;                 % Frequency [rad/s]
a = 5;                      % Gyroscope Bias Factor
b = 2;                      % Gyroscope Bias Factor

% Simulation
n = sigma_n*randn(N,1);     % Noise
r = 100*sin(w*t)';          % 
g = a*r + b*ones(N,1) + n;  % Gyroscope Measurements

% Least Squares
H = [r ones(N,1)];          % 
x = pinv(H)*g;              % State Estimate

fprintf('a) 10 Sample; Individual Least Squares Estimate: [%0.3g %0.3g]\n', x);

% Monte Carlo
Nmc = 1000;
for i = 1:Nmc
    % Simulation
    n = sigma_n*randn(N,1);     % Noise
    g = a*r + b*ones(N,1) + n;  % Gyroscope Measurements

    % Least Squares
    H = [r ones(N,1)];          % 
    x(:,i) = pinv(H)*g;         % State Estimate
end
sigma = std(x,0,2);         % Sample Standard Deviation
x_ = mean(x,2);             % Sample Mean

fprintf('b) 10 Sample; 1000 MC Least Squares Sigma: [%0.3g %0.3g]\n', sigma);
fprintf('b) 10 Sample; 1000 MC Least Squares Estimate: [%0.3g %0.3g]\n\n', x_);

N = 1000;                   % Number of Samples
t = linspace(0,tf,N);       % Time [s]
f = 1/N;                    % Frequency [Hz]
w = 2*pi*f;                 % Frequency [rad/s]

% Simulation
n = sigma_n*randn(N,1);     % Noise
r = 100*sin(w*t)';          % 
g = a*r + b*ones(N,1) + n;  % Gyroscope Measurements

% Least Squares
H = [r ones(N,1)];          % 
x = pinv(H)*g;              % State Estimate

fprintf('c) 1000 Sample; Individual Least Squares Estimate: [%0.3g %0.3g]\n', ...
    x);

% Monte Carlo
Nmc = 1000;
for i = 1:Nmc
    % Simulation
    n = sigma_n*randn(N,1);     % Noise
    g = a*r + b*ones(N,1) + n;  % Gyroscope Measurements

    % Least Squares
    H = [r ones(N,1)];          % 
    x(:,i) = pinv(H)*g;         % State Estimate
end
sigma = std(x,0,2);         % Sample Standard Deviation
x_ = mean(x,2);             % Sample Mean

fprintf('d) 1000 Sample; 1000 MC Least Squares Sigma: [%0.3g %0.3g]\n', sigma);
fprintf('d) 1000 Sample; 1000 MC Least Squares Estimate: [%0.3g %0.3g]\n', x_);

%% PART 4
clear; close all; clc;

N = 1000;                       % Number of Samples

numd = 0.25*[1 -0.8];           % Discrete TF Numerator
dend = [1 -1.9 0.95];           % Discrete TF Denominator
u = randn(N,1);                 % Input Noise
y = dlsim(numd,dend,u);         % Output Simulation
sigma = 0.01;                   % Output Noise Standard Deviation
Y = y + sigma*randn(N,1);       % Output w/ Noise

H = 1;