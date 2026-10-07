function [ stat ] = mc_CovariateCorrectionFaster( Y, X, raw, cols,t)
%MC_COVARIATECORRECTION Correction for a series of covariates using
%multiple regression
%
%   FORMAT [residuals] = mc_CovariateCorrection( Y, X, raw, tvalcalc)
%       Y   -   nExamples x nFeatures matrix of observations
%       X   -   nExamples x nPredictors design matrix
%       raw -   Use this to disable some of the automatic features of mc_CovariateCorrection
%                       0 - Automatically mean center X column-wise and then add intercept
%                       1 - Automatically mean center X, do not add intercept
%                               NOTE - This will mean center columns 2:end, and leave alone column 1 assuming it is intercept
%                       2 - Automatically add intercept, do not mean center
%                       3 - Do not add intercept or mean center
%       cols-   columns to retain betas/t scores from
%       t   -   1 to calculate t scores
%   RESULTS
%       A single struct will be returned. It can have the following fields
%               b       -       nPredicted x nFeatures matrix of betas values from regression
%               t       -       t values corresponding to each b

%
%
% This program will assume the same design matrix for all of your features,
% thus enabling much more rapid computation of the residuals.

X = single(X);
Y = single(Y);

%if(~exist('raw','var') )
%    raw=0;
%end

% Mean center all your covariates
if(raw==0) % if in full helper mode
    X = mc_SweepMean(X); % mean center the whole design matrix
elseif (raw==1) %
    X(:,2:end) = mc_SweepMean(X(:,2:end)); % mean center, but leave first column (intercept)
end

% Prepend a constant to the predictor matrix
if(raw==0 | raw==2)
    X = horzcat(ones(size(X,1),1),X);
end

nFeat = size(Y,2);
nSub = size(Y,1);
nPred = size(X,2);

b =   pinv(X)*Y;
stat.b = b(cols,:);

if (t==1)
    res = Y - X*b;

    C = pinv(X'*X); % X'*X is the covariance matrix of X. The inverse is valuable for denominator or SE(beta hat)

    %xvar_inv = diag(C);
    %xvar_inv = xvar_inv(cols);
    %xvar_inv = repmat(xvar_inv,1,nFeat);

    sse = sum(res.^2,1) ./ (nSub - nPred);
    %sse = repmat(sse,numel(cols),1);

    %bSE = sqrt(xvar_inv .* sse);
    bSE = sqrt(diag(C)*sse);
    
    t = b ./ bSE;
    stat.t = t(cols,:);
    %stat.res = res;
end
