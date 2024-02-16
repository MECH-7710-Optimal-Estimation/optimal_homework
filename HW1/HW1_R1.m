%% QUESTION 1
clear; close all; clc;

N = 6;
a = (1:6)';
b = (4:9)';
c = [1 1 3 3 3 5]';
d1 = a;
d2 = c;
d = [d1 d2];

% a)
[pmf_a, shift_a, mu_a, sigma_a] = simDice(a, N);
% b)
[pmf_b, shift_b, mu_b, sigma_b] = simDice(b, N);
% c)
[pmf_c, shift_c, mu_c, sigma_c] = simDice(c, N);
% d)
[pmf_d, shift_d, mu_d, sigma_d] = sim2Dice(d1, d2, N/2);

figure();
hold('on');
title('Discrete PDFs vs. Continuous Gaussian PDF');
tiledlayout(2,2);
nexttile();
hold('on');
title('1a)');
stem(shift_a, pmf_a);
plot(shift_a, normpdf(shift_a, mu_a, sigma_a));
xlabel('Sum of Dice Rolls');
ylabel('Probability');
legend('Discrete PDF', 'Continuous Gaussian PDF');

nexttile();
hold('on');
title('1b)');
stem(shift_b, pmf_b);
plot(shift_b, normpdf(shift_b, mu_b, sigma_b));
xlabel('Sum of Dice Rolls');
ylabel('Probability');
legend('Discrete PDF', 'Continuous Gaussian PDF');

nexttile();
hold('on');
title('1c)');
stem(shift_c, pmf_c);
plot(shift_c, normpdf(shift_c, mu_c, sigma_c));
xlabel('Sum of Dice Rolls');
ylabel('Probability');
legend('Discrete PDF', 'Continuous Gaussian PDF');

nexttile();
hold('on');
title('1d)');
stem(shift_d, pmf_d);
plot(shift_d, normpdf(shift_d, mu_d, sigma_d));
xlabel('Sum of Dice Rolls');
ylabel('Probability');
legend('Discrete PDF', 'Continuous Gaussian PDF');

fprintf('1a) Mean & Standard Deviation: [%0.3g %0.3g]\n', mu_a, sigma_a)
fprintf('1a) Sum of the PDF: %0.3g\n', sum(pmf_a))
fprintf('1b) Mean & Standard Deviation: [%0.3g %0.3g]\n', mu_b, sigma_b)
fprintf('1b) Sum of the PDF: %0.3g\n', sum(pmf_b))
fprintf('1c) Mean & Standard Deviation: [%0.3g %0.3g]\n', mu_c, sigma_c)
fprintf('1c) Sum of the PDF: %0.3g\n', sum(pmf_c))
fprintf('1d) Mean & Standard Deviation: [%0.3g %0.3g]\n', mu_d, sigma_d)
fprintf('1d) Sum of the PDF: %0.3g\n\n', sum(pmf_d))
%% QUESTION 2
clear;

% a)
x1 = 1:6;
x2 = 1:6;
fx1 = simDice(x1', 1)';
fx2 = simDice(x2', 1)';
fx1x2 = fx1*fx2';
fprintf('2a) fx1x2:\n');
fprintf('\t %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g\n', fx1x2);
Ex1 = x1*fx1x2*ones(1,length(x1))'; % E[x1]
fprintf('\n2a) E[x1]: %0.3g\n', Ex1);
Ex2 = x2*fx1x2*ones(1,length(x2))'; % E[x2]
fprintf('2a) E[x2]: %0.3g\n', Ex2);
Ex1_Ex1 = round((x1 - Ex1)*fx1);    % E[(x1-E[x1])]
fprintf('2a) E[(x1 - E[x1])]: %0.3g\n', Ex1_Ex1);
Ex12 = (x1.^2)*fx1;                 % E[x1^2]
fprintf('2a) E[x1^2]: %0.3g\n', Ex12);
Px1 = ((x1 - Ex1).^2)*fx1;          % E[(x1 - E[x1])^2]
fprintf('2a) E[(x1 - E[x1])^2]: %0.3g\n', Px1);
Px1x2 = ((x1 - Ex1)*(x2 - Ex2)')*fx1x2; % E[(x1 - E[x1])(x2 - E[x2])]
fprintf('2a) Px1x2:\n')
fprintf('\t %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g\n', Px1x2');

% b)
Px1x2 = ((x1 - Ex1)*(x2 - Ex2)')*fx1x2;
fprintf('\n2b) Px1x2:\n')
fprintf('\t %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g\n', Px1x2');

% c)
v1 = x1;
fv1 = fx1;
[fv2, v2] = simDice(x1', 2);
fv2 = fv2';
fv1v2 = fv1*fv2';
fprintf('\n2c) fv1v2:\n');
fprintf(['\t %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g ' ...
    '%0.3g %0.3g %0.3g\n'], fv1v2');

% d)
Ev1 = v1*fv1;
fprintf('\n2d) E[v1]: %0.3g\n', Ev1);
v1RMS = sqrt((v1.^2)*fv1);
fprintf('2d) RMS(v1): %0.3g\n', v1RMS);
Pv1 = ((v1 - Ev1).^2) * fv1;
fprintf('2d) E[(v1 - E[v1])^2): %0.3g\n', Pv1);

% e)
Ev2 = v2*fv2;
fprintf('2e) E[v2]: %0.3g\n', Ev2);
v2RMS = sqrt((v2.^2)*fv2);
fprintf('2e) RMS(v2): %0.3g\n', v2RMS);
Pv2 = ((v2 - Ev2).^2) * fv2;
fprintf('2e) E[(v2 - E[v2])^2): %0.3g\n', Pv2);

% f)
Pv1v2 = ((v1 - Ev1)*(v2(1:6) - Ev2)')*fv1*fv2';
fprintf('2f) Pv1v2:\n');
fprintf(['\t %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g %0.3g ' ...
    '%0.3g %0.3g %0.3g\n'], Pv1v2');

%% QUESTION 4
clear;

Vo = [-2.5 -1.5 -0.5 0.5 1.5 2.5];
fx = groupcounts(Vo')./length(Vo);

figure();
hold('on');
title('4a) Vo PDF')
stem(fx);
xlabel('Value');
ylabel('Probability');

muVo = mean(Vo);
fprintf('\n4b) The mean of Vo: %0.3g\n', muVo);
sigmaVo = std(Vo);
fprintf('4b) The standard deviation of Vo: %0.3g\n', sigmaVo);
varVo = var(Vo);
fprintf('4b) The variance of Vo: %0.3g\n\n', varVo);

%% QUESTION 5
figure(); 
hold('on');
title('5) Probability Function b/w 0 and 2');
subtitle('Probability Equals 0 for All Other Values');
fplot(@(x) x/2, [0 2], '-o');
xline(4/3);
xlabel('Value');
ylabel('Probability');
legend('PDF', 'Mean');

%% QUESTION 6
clear;

Px = [2 1;
      1 4];
[V, D] = eigs(Px);
fprintf('6a) Eigenvalues of Px: [%0.3g %0.3g]\n', diag(D));

c = [0.25, 1, 1.5];
t = linspace(0, 2 * pi);
figure();
hold("on");
title("6c) Likelihood Ellipses for Varying Probabilities")
for k = 1:length(c)
    a = (V * sqrt(c(k)*D)) * [cos(t); sin(t)];
    plot(a(1, :), a(2, :));
end
xlabel("X");
ylabel("Y");
legend('c = 0.25', 'c = 1', 'c = 1.5');

fx = @(c) ((2*pi)^(size(Px,1)/2) * det(Px)^(1/2))^(-1) .* exp(-1/2.*(c.^2));
probs = fx(c);
fprintf('6d) The probability for c = 0.25: %0.3g\n', probs(1));
fprintf('6d) The probability for c = 1: %0.3g\n', probs(2));
fprintf('6d) The probability for c = 1.5: %0.3g\n\n', probs(3));

%% QUESTION 7
clear;
sigmax = 2.0;
varx = sigmax^2;
muy = 2*varx;
vary = 4*3*varx^2 - muy^2;
sigmay = sqrt(vary);

figure();
hold("on");
title('7b) Comparison of PDFs');
plot(-3*sigmay:0.1:3*sigmay, normpdf(-3*sigmay:0.1:3*sigmay, 0, sigmax), '-o');
plot(-3*sigmay:0.1:3*sigmay, normpdf(-3*sigmay:0.1:3*sigmay, muy, sigmay), '-*');
xlabel('X-Value');
ylabel('Probability');
legend('x', 'y = 2x^2');

%% FUNCTIONS
function [pmf, shift, mu, sigma] = simDice(die, N)
    pmf = genPMF(die, N);
    shift = linspace(min(die*N), max(die*N), length(pmf));
    mu = sum(shift.*pmf);
    sigma = sqrt(sum(((shift - mu).^2).*pmf));
end

function [pmf, shift, mu, sigma] = sim2Dice(d1, d2, N)
    pmf = genPMF(d1, N);
    [~, fx] = genPMF(d2, N);
    for i = 1:N
        pmf = conv(pmf, fx);
    end
    pmf(pmf == 0) = [];
    maxs = d1*N + d2*N;
    shift = linspace(min(d1*N + d2*N), max(d1*N + d2*N), length(pmf));
    mu = sum(shift.*pmf);
    sigma = sqrt(sum(((shift - mu).^2).*pmf));
end

function [pmf, fx] = genPMF(die, N)
    [probs, vals] = groupcounts(die);
    fx = zeros(length(die),1);
    fx(vals) = probs./length(die);
    pmf = fx';
    for i = 1:N-1
        pmf = conv(pmf, fx');
    end
    pmf(pmf == 0) = [];
end