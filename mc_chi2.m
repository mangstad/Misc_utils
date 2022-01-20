function [chi2,p_chi2,p_fisher] = mc_chi2(obs_cell, obs_all, cellsize, N)

obs_rest = obs_all - obs_cell;

exp_cell = obs_all * (cellsize/N);
exp_rest = obs_all * ((N-cellsize)/N);

chi2 = ((obs_cell - exp_cell)^2)/exp_cell + ((obs_rest - exp_rest)^2)/exp_rest;
p_chi2 = 1-chi2cdf(chi2,1);

%a = obs_cell;
%b = cellsize - obs_cell;
%c = obs_all-obs_cell;
%d = (N- cellsize) - c;

%p_fisher = (factorial(a+b)*factorial(c+d)*factorial(a+c)*factorial(b+d))/(factorial(a)*factorial(b)*factorial(c)*factorial(d)*factorial(N);

%obs_cell = 50;
%obs_all = 250;
%cellsize = 1000;
%N = 2000;

p_fisher = 1-hygecdf(obs_cell, N, obs_all, cellsize);
