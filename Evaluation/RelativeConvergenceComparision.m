clear
clc

S0 = 50;       % Initial stock price
K = 50;        % Strike price
r = 0.1;       % Risk-free rate
T = 5/12;          % iteration to maturity
sigma = 0.4;    % Volatility
Smax = 100;     % Maximum stock price

iteration = zeros(6, 1);
lu_val = zeros(6, 1);
ja_val = zeros(6, 1);
gs_val = zeros(6, 1);
sor_val = zeros(6, 1);
mul_val = zeros(6, 1);

count = 1;
for i = 50:10:300
    [lu_val(count), ~] = LUSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [ja_val(count), ~] = JacobiSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [gs_val(count), ~] = GaussSeidelSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [sor_val(count), ~] = SORSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [mul_val(count), ~] = MultigridSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    iteration(count) = count;
    count = count + 1;
end

lu_rate = (lu_val(2:length(lu_val)) - lu_val(1:length(lu_val) - 1)) ./ lu_val(1:length(lu_val) - 1);
ja_rate = (ja_val(2:length(ja_val)) - ja_val(1:length(ja_val) - 1)) ./ ja_val(1:length(ja_val) - 1);
gs_rate = (gs_val(2:length(gs_val)) - gs_val(1:length(gs_val) - 1)) ./ gs_val(1:length(gs_val) - 1);
sor_rate = (sor_val(2:length(sor_val)) - sor_val(1:length(sor_val) - 1)) ./ sor_val(1:length(sor_val) - 1);
mul_rate = (mul_val(2:length(mul_val)) - mul_val(1:length(mul_val) - 1)) ./ mul_val(1:length(mul_val) - 1);

figure;
plot(iteration(2:length(iteration)), lu_rate, 'ro-', 'LineWidth', 1); % Red line
hold on;
plot(iteration(2:length(iteration)), ja_rate, 'go-', 'LineWidth', 1); % Green line
plot(iteration(2:length(iteration)), gs_rate, 'bo-', 'LineWidth', 1); % Blue line
plot(iteration(2:length(iteration)), sor_rate, 'mo-', 'LineWidth', 1); % Magenta line
plot(iteration(2:length(iteration)), mul_rate, 'co-', 'LineWidth', 1); % Cyan line
xlabel('iteration');
ylabel('Rate');
title('Relative convergence comparison');
legend('lu\_iteration', 'ja\_iteration', 'gs\_iteration', 'sor\_iteration', 'mul\_iteration');
hold off;
