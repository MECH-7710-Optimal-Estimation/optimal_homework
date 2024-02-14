clear; close all; clc;

N = 6;
Nd = 3;

a = [1 2 3 4 5 6]';
b = [4 5 6 7 8 9]';
c = [1 1 3 3 3 5]';
d1 = [1 2 3 4 5 6]';
d2 = [1 1 3 3 3 5]';

[pdf_a, shift_a] = dice_sim(a, N);
[mu_a, sigma_a] = dice_stats(a, N);

[pdf_b, shift_b] = dice_sim(b, N);
[mu_b, sigma_b] = dice_stats(b, N);

[pdf_c, shift_c] = dice_sim(c, N);
[mu_c, sigma_c] = dice_stats(c, N);

[pdf_d1, shift_d1] = dice_sim(d1, N);
[mu_d1, sigma_d1] = dice_stats(d1, N);
[pdf_d2, shift_d2] = dice_sim(d2, N);
[mu_d2, sigma_d2] = dice_stats(d2, N);
pdf_d = conv(pdf_d1, pdf_d2);
big_shift = vertcat(shift_d1', shift_d2');
shift_d = linspace(min(big_shift), max(big_shift), length(pdf_d));

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

figure();
stem(shift_d, pdf_d);
fprintf('The Sum of the PDF for d is : %0.3g\n', sum(pdf_d));

%% QUESTION 4
clear;
Vo = [-2.5 -1.5 -0.5 0.5 1.5 2.5];

mu = mean(Vo);
sig = std(Vo);
plot((mu - 3*sig):0.001:(mu + 3*sig), ...
    normpdf((mu - 3*sig):0.001:(mu + 3*sig), mu, sig));

idx = 1;
for r = 0:0.1:1
    Vn = zeros(length(Vo),1);
    for k = 1:length(Vo)
        Vn(k+1) = (1-r)*Vn(k) + r*Vo(k);
    end
    m(idx) = mean(Vn);
    s(idx) = std(Vn);
    idx = idx + 1;
end

figure();
hold("on");
plot(0:0.1:1, m);
plot(0:0.1:1, s);

%% QUESTION 5
clear;

figure();
fplot(@(x) x/2, [0 2])

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