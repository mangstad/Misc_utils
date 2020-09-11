function results = mc_bbs(featuremat,pheno,nuisance,folds,NumComp,NumPerms,varargin)

p = inputParser;
% addRequired(p,'featuremat');
% addRequired(p,'pheno');
% addRequired(p,'nuisance');
% addRequired(p,'folds');
% addRequired(p,'NumComp');
% addRequired(p,'NumPerms');
addParameter(p,'Scores',[]);
addParameter(p,'Components',[]);
addParameter(p,'LOSOPheno',0);
parse(p,varargin{:});

good = ~any(isnan(pheno),2) & ~any(isnan(nuisance),2);
featuremat = featuremat(good,:);
pheno = pheno(good,:);
nuisance = nuisance(good,:);
folds = folds(good);

n = size(featuremat,1);

nFold = numel(unique(folds));

pcadone = 0;
if (~isempty(p.Results.Scores))
    Aa = p.Results.Scores;    
    if (size(Aa{1},1)>n)
        for iFold = 1:nFold
            Aa{iFold} = Aa{iFold}(good,:);
        end
    end
    pcadone = 1;
    results.Aa = Aa;
end
if (~isempty(p.Results.Components))
    components = p.Results.Components;
    results.components = components;
end
LOSOPheno = p.Results.LOSOPheno;

% if (nargin > 6) 
%     Aa = varargin{1};
%     components = varargin{2};
%     pcadone = 1;
%     if (size(Aa{1},1)>n)
%         for iFold = 1:nFold
%             Aa{iFold} = Aa{iFold}(good,:);
%         end
%     end
% end

good = good(good);

if (pcadone==0) 
    %clear Aa components;
    for iFold = 1:nFold
        fprintf(1,'.');
        %find train and test data for this fold
        test_idx = folds==iFold & good;
        train_idx = ~test_idx & good;
        
        %reduce the training data
        coeff = pca(featuremat(train_idx,:));
        components{iFold} = coeff';

        %mean center train, and mean center test with train means
        mu = mean(featuremat(train_idx,:));
        x = bsxfun(@minus,featuremat,mu);

        %calculate expressions for each subject for train and test
        Aa{iFold} = (pinv(components{iFold}')*x')';
    end
    results.Aa = Aa;
    results.components = components;
end

pheno_predict = zeros(n,size(pheno,2));
pheno_residualized = zeros(n,size(pheno,2));

for iFold = 1:nFold
    %find train and test data for this fold
    test_idx = folds==iFold & good;
    train_idx = ~test_idx & good;
    
    n_trainf = sum(train_idx);
    n_testf = sum(test_idx);

    k = NumComp;

    Abig = Aa{iFold}(train_idx,1:k);
    Abig_test = Aa{iFold}(test_idx,1:k);
    %predicting phenotype
    for iPheno = 1:size(pheno,2)
        X = [ones(n_trainf,1) Abig nuisance(train_idx,:)];
        b = pinv(X'*X)*X'*pheno(train_idx,iPheno);

        pheno_predict(test_idx,iPheno) = Abig_test*b(2:(k+1));
        pheno_residualized(test_idx,iPheno) = pheno(test_idx,iPheno) - [ones(n_testf,1) nuisance(test_idx,:)]*b([1 (k+2):end]);
    end
        
end

%check correlation between actual and predicted phenotypes
fold_corr = zeros(size(pheno,2),nFold);
if (LOSOPheno)
    fold_corr = zeros(1,nFold);
end

for iFold = 1:nFold
    test_idx = folds==iFold & good;
    if (LOSOPheno)
        fold_corr(:,iFold) = diag(corr(pheno_predict(test_idx,iFold),pheno_residualized(test_idx,iFold)));    
    else
        fold_corr(:,iFold) = diag(corr(pheno_predict(test_idx,:),pheno_residualized(test_idx,:)));
    end
end

mean_corr = mc_FisherZ(mean(mc_FisherZ(fold_corr),2),1);
std_corr = mc_FisherZ(std(mc_FisherZ(fold_corr),[],2),1);
ts_corr = tinv(0.975,nFold-1);
ci_corr = ts_corr.*std_corr./sqrt(nFold);

results.pheno_predict = pheno_predict;
results.pheno_residualized = pheno_residualized;
results.fold_corr = fold_corr;
results.mean_corr = mean_corr;
results.std_corr = std_corr;
results.ts_corr = ts_corr;
results.ci_corr = ci_corr;


shuf_idx = zeros(n,NumPerms);
for iPerm = 1:NumPerms
    %shuffle subjects within site
    for iFold = 1:nFold
        currsite = find(folds==iFold);
        shuf_idx(currsite,iPerm) = randsample(currsite,numel(currsite));
    end
end

%check correlation between actual and predicted phenotypes
fold_corr_perm = zeros(size(pheno,2),nFold,NumPerms);
mean_corr_perm = zeros(size(pheno,2),NumPerms);
if (LOSOPheno==1) 
    fold_corr_perm = zeros(1,nFold,NumPerms);
    mean_corr_perm = zeros(1,NumPerms);
end

P = size(pheno,2);

parfor iPerm = 1:NumPerms
    fprintf(1,'%d of %d\n',iPerm,NumPerms);
    pheno_predict_perm = zeros(n,size(pheno,2));
    pheno_residualized_perm = zeros(n,size(pheno,2));

    for iFold = 1:nFold
        %find train and test data for this fold
        test_idx = folds==iFold & good;
        train_idx = ~test_idx & good;

        n_trainf = sum(train_idx);
        n_testf = sum(test_idx);

        k = NumComp;

        Abig = Aa{iFold}(train_idx,1:k);
        Abig_test = Aa{iFold}(test_idx,1:k);

        X = [ones(n,1) nuisance];
        b = pinv(X(train_idx,:)'*X(train_idx,:))*X(train_idx,:)'*pheno(train_idx,:);

        %get residuals and permute
        res = pheno - X*b;
        res = res(shuf_idx(:,iPerm),:);

        %add true nuisance effects back to permuted residuals
        shufpheno = X*b + res;

        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            X = [ones(n_trainf,1) Abig nuisance(train_idx,:)];
            b = pinv(X'*X)*X'*shufpheno(train_idx,iPheno);

            pheno_predict_perm(test_idx,iPheno) = Abig_test*b(2:(k+1));
            pheno_residualized_perm(test_idx,iPheno) = shufpheno(test_idx,iPheno) - [ones(n_testf,1) nuisance(test_idx,:)]*b([1 (k+2):end]);
        end     
    end
    
    tmp = zeros(P,nFold);
    for iFold = 1:nFold
        test_idx = folds==iFold & good;    
        if (LOSOPheno)
            tmp(:,iFold) = diag(corr(pheno_predict_perm(test_idx,iFold),pheno_residualized_perm(test_idx,iFold)));    
        else
            tmp(:,iFold) = diag(corr(pheno_predict_perm(test_idx,:),pheno_residualized_perm(test_idx,:)));
        end
    end

    fold_corr_perm(:,:,iPerm) = tmp;
    mean_corr_perm(:,iPerm) = squeeze(mc_FisherZ(mean(mc_FisherZ(tmp),2),1));
end

results.fold_corr_perm = fold_corr_perm;
results.mean_corr_perm = mean_corr_perm;

%calculate perm p values
results.perm_p = (1+sum(bsxfun(@gt,results.mean_corr_perm,results.mean_corr),2))./NumPerms;
