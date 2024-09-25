function [r,ci] = mc_corr_ci(x,y,varargin)
    conf = 0.95;
    if (nargin>2)
        conf = varargin{1};
    end
    alpha = 1-conf;
    
    n = size(x,1);
    r = corr(x,y);
    z = mc_FisherZ(r);
    z1a = icdf('norm',1-(alpha/2),0,1);
    ci(1) = mc_FisherZ(z - (z1a/sqrt(n-3)),1);
    ci(2) = mc_FisherZ(z + (z1a/sqrt(n-3)),1);
    