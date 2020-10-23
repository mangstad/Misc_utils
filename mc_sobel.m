function [t,p,coefs] = mc_sobel(X,M,Y)

%X = zscore(X);
%M = zscore(M);
%Y = zscore(Y);

m1 = fitlm(X,Y);
m2 = fitlm(X,M);
m3 = fitlm([X M],Y);

cp = m3.Coefficients.Estimate(2);
scp = m3.Coefficients.SE(2);

c = m1.Coefficients.Estimate(2);
sc = m1.Coefficients.SE(2);

a = m2.Coefficients.Estimate(2);
sa = m2.Coefficients.SE(2);

b = m3.Coefficients.Estimate(3);
sb = m3.Coefficients.SE(3);

se = sqrt(a^2*sb^2 + b^2*sa^2);

t = (a*b)/se;

p = 2*(normcdf(-abs(t),0,1));

coefs = [a b c cp];
