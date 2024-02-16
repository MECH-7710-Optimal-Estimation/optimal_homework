clear; close all; clc;

N = 6;
Nd = 3;

a = [1 2 3 4 5 6]';
b = [4 5 6 7 8 9]';
c = [1 1 3 3 3 5]';
d1 = [1 2 3 4 5 6]';
d2 = [1 1 3 3 3 5]';

[pdf_a, shift_a] = dice_sim(a, N);
% [mu_a, sigma_a] = dice_stats(a, N);
mu_a = sum(shift_a*pdf_a);
sigma_a = sqrt((shift_a.^2) * pdf_a);

[pdf_b, shift_b] = dice_sim(b, N);
[mu_b, sigma_b] = dice_stats(b, N);

[pdf_c, shift_c] = dice_sim(c, N);
[mu_c, sigma_c] = dice_stats(c, N);

% [pdf_d1, shift_d1] = dice_sim(d1, N);
% [mu_d1, sigma_d1] = dice_stats(d1, N);
% [pdf_d2, shift_d2] = dice_sim(d2, N);
% % [mu_d2, sigma_d2] = dice_stats(d2, N);
% pdf_d = conv(pdf_d1, pdf_d2);
% % mu_d = mu_d1 + mu_d2;
% % sigma_d = sqrt(sigma_d1^2 + sigma_d2^2);
% big_shift = vertcat(shift_d1', shift_d2');
% shift_d = linspace(min(big_shift), max(big_shift), length(pdf_d));

pmf = groupcounts(d1)./length(d1)

pdf = pmf;
for i = 1:2
    pdf = conv(pdf, pmf)
end

pmf = groupcounts(d2)./length(d2)
pmf = [pmf(1) 0 pmf(2) 0 pmf(3) 0]';
for i = 1:3
    pdf = conv(pdf, pmf)
end

x = (6:33)';
pdf = pdf(1:end-3);
mu = sum(x.*pdf)
sig = sqrt(sum(((x - mu).^2).*pdf))

% sig = sqrt(sum((((6:36) - mu).^2) .* pdf));

figure();
hold on
stem(x, pdf);
plot(x, normpdf(x, mu, sig))



figure();
hold("on");
stem(shift_a, pdf_a);
plot(shift_a, normpdf(shift_a, mu_a, sigma_a));

fprintf('The Sum of the PDF for a is : %0.3g\n', sum(pdf_a));

figure();
hold("on");
stem(shift_b, pdf_b);
plot(shift_b, normpdf(shift_b, mu_b, sigma_b));

fprintf('The Sum of the PDF for b is : %0.3g\n', sum(pdf_b));

figure();
hold("on");
stem(shift_c, pdf_c);
plot(shift_c, normpdf(shift_c, mu_c, sigma_c));

fprintf('The Sum of the PDF for c is : %0.3g\n', sum(pdf_c));

% figure();
% hold("on")
% stem(shift_d, pdf_d);
% fprintf('The Sum of the PDF for d is : %0.3g\n', sum(pdf_d));
% plot(shift_d, normpdf(shift_d, mu_d, sigma_d));

%% QUESTION 2
clear;

% a)
x1 = 1:6;
x2 = 1:6;
fx1 = dice_sim(x1', 1);
fx2 = dice_sim(x2', 1);
fx1x2 = fx1*fx2';
Ex1 = x1*fx1x2*ones(1,length(x1))';
Ex2 = x2*fx1x2*ones(1,length(x2))';
Ex1_x1_ = round((x1 - Ex1)*fx1);
Ex12 = (x1.^2)*fx1;
Ex1_x1_2 = ((x1 - Ex1).^2)*fx1;
Ex1_x1_x2_x2_ = ((x1 - Ex1)*(x2 - Ex2)')*fx1x2;

% b)
Px1x2 = ((x1 - Ex1)*(x2 - Ex2)')*fx1x2;

% c)
v1 = x1;
fv1 = fx1;
[fv2, v2] = dice_sim(x1', 2);
fv1v2 = fv1*fv2';

% d)
Ev1 = v1*fv1;
Pv1 = ((v1 - Ev1).^2) * fv1;

% e)
Ev2 = v2*fv2;
Pv2 = ((v2 - Ev2).^2) * fv2;
Pv1v2 = ((v1 - Ev1)*(v2(1:6) - Ev2)')*fv1*fv2';

%% QUESTION 4
clear;
Vo = [-2.5 -1.5 -0.5 0.5 1.5 2.5]' + 2*ones(6,1);
pmf = grouptransform(Vo, Vo, @numel)./length(Vo);

% figure();
% stem(min(Vo):max(Vo), pmf);

muVo = mean(Vo);
varVo = var(Vo);

r = 0.25;
figure();
hold("on");
idx = 1;
for r = 0.1:1e-1:1
    Vn = 3*ones(100,1);
    for i = 1:1000
        for k = 1:length(Vn)-1
            Vn(k + 1) = (1-r)*Vn(k) + r*Vo(randi(length(Vo)));
        end
        m(idx) = mean(Vn(round(3*length(Vn)/4):end));
    end
    idx = idx +1;
    plot(Vn)
end

mean(m)
% plot((mu - 3*sig):0.001:(mu + 3*sig), ...
%     normpdf((mu - 3*sig):0.001:(mu + 3*sig), mu, sig));

% idx = 1;
% for r = 0:0.1:1
%     Vn = zeros(length(Vo),1);
%     for k = 1:length(Vo)
%         Vn(k+1) = (1-r)*Vn(k) + r*Vo(k);
%     end
%     m(idx) = mean(Vn);
%     s(idx) = std(Vn);
%     idx = idx + 1;
% end

% figure();
% hold("on");
% plot(0:0.1:1, m);
% plot(0:0.1:1, s);

%% QUESTION 5
clear;

figure();
fplot(@(x) x/2, [0 2])
xline(4/3)

%% QUESTION 6
clear;

Px = [2 1;
      1 4];
[V, D] = eigs(Px);
lambda = diag(D);
Pa1 = (1/sqrt(lambda(1))) * (V(:,1)/norm(V(:,1)));
Pa2 = (1/sqrt(lambda(2))) * (V(:,2)/norm(V(:,2)));

c = [0.25, 1, 1.5];
t = linspace(0, 2 * pi);
figure();
hold("on");
title("Likelihood Ellipses for Varying Probabilities")
for k = 1:length(c)
    a = (V * sqrt(c(k)*D)) * [cos(t); sin(t)];
    plot(a(1, :), a(2, :));
end
xlabel("X");
ylabel("Y");
legend('c = 0.25', 'c = 1', 'c = 1.5');

fx = @(c) ((2*pi)^(size(Px,1)/2) * det(Px)^(1/2))^(-1) .* exp(-1/2.*(c.^2));
probs = fx(c);
fprintf('The probability for c = 0.25: %0.3g\n', probs(1));
fprintf('The probability for c = 1: %0.3g\n', probs(2));
fprintf('The probability for c = 1.5: %0.3g\n\n', probs(3));

%% QUESTION 7
clear;
sigmax = 2.0;
varx = sigmax^2;
muy = 2*varx;
vary = 4*3*varx^2 - muy^2;
sigmay = sqrt(vary);

figure();
hold("on");
plot(-3*sigmay:0.1:3*sigmay, normpdf(-3*sigmay:0.1:3*sigmay, 0, sigmax));
plot(-3*sigmay:0.1:3*sigmay, normpdf(-3*sigmay:0.1:3*sigmay, muy, sigmay));

fprintf(['The mean is 2*sigmax^2 which is essentially plugging sigma \n' ...
    'into y(x).  Sigmay is larger than sigmax.\n']);
fprintf('y is a normal random variable\n\n');

%% FUNCTIONS
function [pdf, shift_range] = dice_sim(die, N)
    pmf = groupcounts(die)./length(die);
    max_vals = conv(die, N);
    
    pdf = pmf;
    for i = 1:N-1
        pdf = conv(pdf, pmf);
    end
    shift_range = linspace(min(max_vals), max(max_vals), length(pdf));
end

function [mu, sigma] = dice_stats(die, N)
    pmf = grouptransform(die, die, @numel)./length(die);
    [~,c] = unique(die, 'first');
    z = zeros(length(die),1);
    z(c) = 1;
    pmf(z == 0) = 0;
    mu = sum(die.*pmf);
    sigma2 = N*sum(((die - mu).^2).*pmf);
    mu = N*mu;
    sigma = sqrt(sigma2);
end