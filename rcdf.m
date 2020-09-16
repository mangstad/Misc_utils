function p = rcdf(R,N)
t = R.*sqrt((N-2)./(1-R.^2));
p = 2*tcdf(-abs(t),N-2);

