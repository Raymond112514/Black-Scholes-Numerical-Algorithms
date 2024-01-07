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


figure;
plot(iteration, lu_val, 'ro-', 'LineWidth', 1); % Red line
hold on;
plot(iteration, ja_val, 'go-', 'LineWidth', 1); % Green line
plot(iteration, gs_val, 'bo-', 'LineWidth', 1); % Blue line
plot(iteration, sor_val, 'mo-', 'LineWidth', 1); % Magenta line
plot(iteration, mul_val, 'co-', 'LineWidth', 1); % Cyan line
xlabel('iteration');
ylabel('Values');
title('Convergence comparison');
legend('lu\_iteration', 'ja\_iteration', 'gs\_iteration', 'sor\_iteration', 'mul\_iteration');
hold off;