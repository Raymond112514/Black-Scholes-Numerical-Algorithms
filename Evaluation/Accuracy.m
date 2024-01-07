clear
clc

lu_error = zeros(6, 1);
ja_error = zeros(6, 1);
gs_error = zeros(6, 1);
sor_error = zeros(6, 1);
mul_error = zeros(6, 1);

for numIter = 1:5
    S0 = 50;       % Initial stock price
    K = 50;        % Strike price
    r = 0.5 * rand();       % Risk-free rate
    T = rand();          % iteration to maturity
    sigma = 1.2 * rand();    % Volatility
    Smax = 100;     % Maximum stock price
    size = 250;
    
    [lu_val, ~] = LUSolver(S0, K, r, T, sigma, Smax, size, size, false, false);
    [ja_val, ~] = JacobiSolver(S0, K, r, T, sigma, Smax, size, size, false, false);
    [gs_val, ~] = GaussSeidelSolver(S0, K, r, T, sigma, Smax, size, size, false, false);
    [sor_val, ~] = SORSolver(S0, K, r, T, sigma, Smax, size, size, false, false);
    [mul_val, ~] = MultigridSolver(S0, K, r, T, sigma, Smax, size, size, false, false);
    truth = blackScholesPut(S0, K, r, T, sigma);
    
    lu_error(numIter) = (lu_val - truth).^2;
    ja_error(numIter) = (ja_val - truth).^2;
    gs_error(numIter) = (gs_val - truth).^2;
    sor_error(numIter) = (sor_val - truth).^2;
    mul_error(numIter) = (mul_val - truth).^2;
   
end





