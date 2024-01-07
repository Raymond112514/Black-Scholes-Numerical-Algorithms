clear
clc

S0 = 50;       % Initial stock price
K = 50;        % Strike price
r = 0.1;       % Risk-free rate
T = 5/12;          % Time to maturity
sigma = 0.4;    % Volatility
Smax = 100;     % Maximum stock price

%4.0760

time = zeros(6, 1);
lu_time = zeros(6, 1);
ja_time = zeros(6, 1);
gs_time = zeros(6, 1);
sor_time = zeros(6, 1);
mul_time = zeros(6, 1);

count = 1;
for i = 50:10:200
    [~, lu_time(count)] = LUSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [~, ja_time(count)] = JacobiSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [~, gs_time(count)] = GaussSeidelSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [~, sor_time(count)] = SORSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    [~, mul_time(count)] = MultigridSolver(S0, K, r, T, sigma, Smax, i, i, false, false);
    time(count) = i;
    count = count + 1;
end


figure;
plot(time, lu_time, 'ro-', 'LineWidth', 1); % Red line
hold on;
plot(time, ja_time, 'go-', 'LineWidth', 1); % Green line
plot(time, gs_time, 'bo-', 'LineWidth', 1); % Blue line
plot(time, sor_time, 'mo-', 'LineWidth', 1); % Magenta line
plot(time, mul_time, 'co-', 'LineWidth', 1); % Cyan line
xlabel('Problem size');
ylabel('Time (sec)');
title('Time complexity of different algorithms');
legend('lu\_time', 'ja\_time', 'gs\_time', 'sor\_time', 'mul\_time');
hold off;
