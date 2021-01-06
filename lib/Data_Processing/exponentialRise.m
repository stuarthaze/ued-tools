function Y = exponentialRise(X,A,B,X0)
    % X is an array of x-values
    % A is a vector of amplitudes
    % B is a vector of time constants
    % X0 is zero position
    nexp = length(A);
    Y = zeros(size(X));
    for ii = 1:nexp
        Y = Y + A(ii)*(1-exp(-(X-X0)/B(ii)));
    end
    % set values at x <= 0 to 0
    Y(X <= X0) = 0;    
end