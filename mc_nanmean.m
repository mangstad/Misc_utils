function mean = mc_nanmean(X)
    index = isnan(X);
    Y = X;
    Y(index) = 0;
    n = sum(~index,1);
    n(n==0) = NaN;
    mean = sum(Y,1)./n;
    