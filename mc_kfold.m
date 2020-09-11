function folds = mc_kfold(n,k,varargin)

p = inputParser;
%addRequired(p,'n');
%addRequired(p,'k');
addParameter(p,'Group',[]);
addParameter(p,'Seed',[]);
parse(p,varargin{:});

Group = p.Results.Group;
Seed = p.Results.Seed;

if (~isempty(Seed))
    rng(Seed);
end

if (~isempty(Group))
    %switch things around for group k-fold
    u = unique(Group);
    n_orig = n;
    n = numel(u);
end

if (n==k)
    folds = randsample(1:k,n,0);
else
    %fold_sizes = repmat(floor(n/k),k,1);
    %fold_sizes(1:mod(n,k)) = fold_sizes(1:mod(n,k))+1;
    %build vector of fold indicies
    tmp = [kron([1:k]',ones(k,1));[1:mod(n,k)]'];
    tmp = [repelem(1:k,floor(n/k))';[1:mod(n,k)]'];
    %permute that vector
    folds = tmp(randperm(n));

end

if (~isempty(Group))
    groupfolds = folds;
    folds = zeros(n_orig,1);
    for i = 1:numel(u)
        folds(Group==u(i)) = groupfolds(i);
    end
end
