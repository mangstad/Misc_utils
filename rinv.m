function r = rinv(P,N)
t = tinv(1-(P/2),N-2);
r = sqrt((t.^2)./(t.^2 + (N-2)));

