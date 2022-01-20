function [cfolds,mapping] = mc_correct_folds(folds)
%this function just renumbers folds so that they are 1:nFold
u = unique(folds);
nFold = numel(u);
cfolds = 0*folds;
for i = 1:nFold
    cfolds(folds==u(i)) = i;
end

mapping = [u [1:nFold]'];
