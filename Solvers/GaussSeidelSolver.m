% Black Scholes PDE (for put options) implicit method solver using Gauss Seidel method
% Input:
%   S0: Initial stock price
%   K: Strike price
%   r: Risk-free rate
%   T: Time to maturity
%   sigma: Volatility
%   Smax: Maximum stock price
%   M, N: Number of steps for price and time
%   message: If set to true, final price and time elaspsed is displayed
%   plot: If set to true, the surface for the option price is displayed

function [price, elapsedTime] = GaussSeidelSolver(S0, K, r, T, sigma, Smax, M, N, text, verbose)
    
    % Start the timer
    tic;
    
    % Initialize relevant parameters
    dt = T / N;
    
    % Construct mesh matrix and set up boundary conditions
    mesh = zeros(M + 1, N + 1);
    S = linspace(0, Smax, M + 1);
    veti = 0:M;
    vetj = 0:N;
    mesh(:, N + 1) = max(K - S, 0);
    mesh(1, :) = K * exp(-r * dt * (N - vetj));
    mesh(M + 1, :) = 0;
    
    % Construct the coefficient matrix 
    a = 0.5 * (r * dt * veti - sigma^2 * dt * (veti.^2));
    b = 1 + sigma^2 * dt * (veti.^2) + r * dt;
    c = -0.5 * (r * dt * veti + sigma^2 * dt * (veti.^2));
    coeff = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1), 1);
    
    % Parameters for Gauss Seidel method
    aux = zeros(M-1 ,1);  
    tol = 10.^-9;
    maxIter = 100;

    % Perform Gauss Seidel iteration from N to 1
    for j = N:-1:1
        aux(1) = -a(2) * mesh(1, j);
        aux(end) = -c(end) * mesh(end, j);
        target = mesh(2:M, j+1) + aux;
        for iter = 1:maxIter
            xold = mesh(2:M, j);
            for i = 1:M-1
                sum1 = dot(coeff(i, 1:i - 1), mesh(2:i, j));
                sum2 = dot(coeff(i, i + 1:end), mesh(i+2:M, j));
                mesh(i + 1, j) = (target(i) - sum1 - sum2) / coeff(i, i);
            end
            if norm(xold - mesh(2:M, j)) < tol
                break;
            end
        end
    end
    
    % Perform interpolation to get price
    price = interp1(S, mesh(:, 1), S0);
    elapsedTime = toc;
    
    % Print out output
    if text
        disp(['Number of meshpoints: (', num2str(M), ',', num2str(N), ')']);
        disp(['Option price: ', num2str(price), ' dollars']);
        disp(['Elapsed Time: ', num2str(elapsedTime), ' seconds']);
    end
    
    % Plot the resuling surface
    if verbose
        figure;
        surf(mesh);
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Function Value');
        title('Surface Plot');
    end
end