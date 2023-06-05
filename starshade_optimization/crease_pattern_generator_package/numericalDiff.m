function df = numericalDiff(f, x, scheme)
    n = length(x);
    if n ~= length(f)
        display 'f and x have unequal lengths'
    end
    
    df = nan(n,1);
    
    switch scheme
        case 'forward1'
            df(1:end-1) = (f(2:end) - f(1:end-1))./(x(2:end) - x(1:end-1));
            df(end) = df(end-1);
        case 'backward1'
            df(2:end) = (f(2:end) - f(1:end-1))./(x(2:end) - x(1:end-1));
            df(1) = df(2);
        case 'central2'
            df(2:end-1) = (f(3:end) - f(1:end-2))./(x(3:end) - x(1:end-2));
            df_init = numericalDiff(f(1:2),x(1:2),'forward1');
            df_end = numericalDiff(f(end-1:end),x(end-1:end),'backward1');
            df(1) = df_init(1);
            df(end) = df_end(end);
    end
end