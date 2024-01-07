function putOptionPrice = blackScholesPut(S0, K, r, T, sigma)
    % Calculate d1 and d2
    d1 = (log(S0 / K) + (r + (sigma^2) / 2) * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);

    % Use cumulative distribution function (normcdf) to compute N(-d1) and N(-d2)
    N_d1 = 0.5 * (1 + erf(-d1 / sqrt(2)));
    N_d2 = 0.5 * (1 + erf(-d2 / sqrt(2)));

    % Calculate put option price using the Black-Scholes formula
    putOptionPrice = K * exp(-r * T) * N_d2 - S0 * N_d1;
end