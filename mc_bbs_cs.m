function results = mc_bbs_cs(featuremat,pheno,nuisance,folds,NumComp,varargin)

n = size(featuremat,1);

pcadone = 0;
if (nargin > 5) 
    Aa = varargin{1};
    components = varargin{2};
    pcadone = 1;
end

thresh = 1;
if (nargin > 7)
    thresh = varargin{3};
end

good = ~any(isnan(pheno),2) & ~any(isnan(nuisance),2);

nFold = numel(unique(folds));
if (pcadone==0) 
clear Aa components;
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
end

results.Aa = Aa;
results.components = components;

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
        mdl = fitlm([Abig nuisance(train_idx,:)],pheno(train_idx,iPheno));
        p = mdl.Coefficients.pValue(2:(k+1));
        idx = p<thresh;
        ks = sum(idx);
        X = [ones(n_trainf,1) Abig(:,idx) nuisance(train_idx,:)];
        b = pinv(X'*X)*X'*pheno(train_idx,iPheno);
    
        pheno_predict(test_idx,iPheno) = Abig_test(:,idx)*b(2:(ks+1));
        pheno_residualized(test_idx,iPheno) = pheno(test_idx,iPheno) - [ones(n_testf,1) nuisance(test_idx,:)]*b([1 (ks+2):end]);
    end
        
end

%check correlation between actual and predicted phenotypes
fold_corr = zeros(size(pheno,2),nFold);
for iFold = 1:nFold
    test_idx = folds==iFold & good;
    fold_corr(:,iFold) = diag(corr(pheno_predict(test_idx,:),pheno_residualized(test_idx,:)));
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
