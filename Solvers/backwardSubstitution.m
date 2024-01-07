function x = backwardSubstitution(U, y)
    
    n = length(y);
    x = zeros(n, 1);

    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:end) * x(i+1:end)) / U(i, i);
    end
end

