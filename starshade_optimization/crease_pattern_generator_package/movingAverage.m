% returns the moving average mu of a series of scalars x or a series of
% column vectors x

function [mu, sigma] = movingAverage(X, span)

    % want span to be odd
    if mod(span, 2) == 0
        span = span - 1;
    end       
    
    N       = max(size(X)); % number of points in the series
    M       = min(size(X)); % dimension of the vector X; expect M < N
    if M > N
        display('movingAverage may not work as expected since M > N')
    end
    
    X       = reshape(X, M, N);
    mu      = nan(M, N);
    sigma   = nan(M, N);
    
    for i = 1:N
        if (2*i - 1) < span
            pencil = 1:(2*i - 1);
        elseif (2*(N - i) + 1) < span
            pencil = (N - 2*(N - i)):N;
        else
            pencil = (i - (span - 1)/2):(i  + (span - 1)/2);
        end
        
        mu(:, i)    = sum(X(:, pencil), 2)/length(pencil);
        sigma(:, i) = sqrt(sum((X(:, pencil) - mu(i)).^2, 2)/length(pencil));        
    end   
end