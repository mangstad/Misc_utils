function results = mc_bbs_logistic(featuremat,pheno,nuisance,folds,NumComp,varargin)

p = inputParser;
% addRequired(p,'featuremat');
% addRequired(p,'pheno');
% addRequired(p,'nuisance');
% addRequired(p,'folds');
% addRequired(p,'NumComp');
addParameter(p,'Scores',[]);
addParameter(p,'Components',[]);
addParameter(p,'LOSOPheno',0);
addParameter(p,'TestPheno',[]);

parse(p,varargin{:});

results.Aa = p.Results.Scores;
results.components = p.Results.Components;
Aa = results.Aa;

n = size(featuremat,1);

nFold = numel(unique(folds));

NumComps = numel(NumComp);

testpheno = p.Results.TestPheno;

pcadone = 1;

LOSOPheno = p.Results.LOSOPheno;

clear p;

good = ~any(isnan(pheno),2) & ~any(isnan(nuisance),2);

nFold = numel(unique(folds));

pheno_predict = zeros(n,size(pheno,2));
pheno_residualized = pheno;
if (~isempty(testpheno))
    pheno_residualized = testpheno;
else
    
end

for iFold = 1:nFold
    %find train and test data for this fold
    test_idx = folds==iFold & good;
    train_idx = ~test_idx & good;
    
    n_trainf = sum(train_idx);
    n_testf = sum(test_idx);

    k = NumComp;

    %if NumComps>1 then do nested 10-fold
    if (NumComps>1)
        nestfold = randsample(10,n_trainf,1);
        nestAa = [];
        for iNest = 1:10
            nestAa{iNest} = Aa{iFold}(train_idx,:);
        end
        nestresults = zeros(size(pheno,2),NumComps);
        for iNest = 1:NumComps
            tempresults = mc_bbs_logistic(featuremat(train_idx,:),pheno(train_idx,:),nuisance(train_idx,:),nestfold,NumComp(iNest),'Scores',nestAa,'LOSOPheno',LOSOPheno);
            nestresults(:,iNest) = tempresults.mean_corr;
        end
        [bestresults,bestcomps] = max(nestresults,[],2);
        bestcomps = NumComp(bestcomps)';
        results.bestcomps(:,iFold) = bestcomps;
        Abig = Aa{iFold}(train_idx,1:max(NumComp));
        Abig_test = Aa{iFold}(test_idx,1:max(NumComp));
        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            kk = bestcomps(iPheno);
            X = [ones(n_trainf,1) Abig(:,1:kk) nuisance(train_idx,:)];
            b = pinv(X'*X)*X'*pheno(train_idx,iPheno);

            pheno_predict(test_idx,iPheno) = Abig_test(:,1:kk)*b(2:(kk+1));
        end
    else
        Abig = Aa{iFold}(train_idx,1:k);
        Abig_test = Aa{iFold}(test_idx,1:k);
        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            %calculate nuisance only model to cleanse train and test
            X = [ones(n_trainf,1) nuisance(train_idx,:)];
            b = pinv(X'*X)*X'*Abig;
            Abig_res = Abig - X*b;
            X_test = [ones(n_testf,1) nuisance(test_idx,:)];
            Abig_test_res = Abig_test - X_test*b;
            
            mdl = fitglm(Abig_res,pheno(train_idx,iPheno),'Distribution','binomial','Link','logit');
            pheno_predict(test_idx,iPheno) = predict(mdl,Abig_test_res);
        end
    end
end

%check correlation between actual and predicted phenotypes
fold_auc = zeros(size(pheno,2),nFold);

if (LOSOPheno)
    fold_auc = zeros(1,nFold);
end

for iFold = 1:nFold
    test_idx = folds==iFold & good;
    train_idx = ~test_idx & good;
    if (LOSOPheno)
        [~,~,~,fold_auc(:,iFold)] = perfcurve(pheno_residualized(test_idx,iFold),pheno_predict(test_idx,iFold),1);   
    else
        for iPheno = 1:size(pheno,2)
            try
            [~,~,~,fold_auc(iPheno,iFold)] = perfcurve(pheno_residualized(test_idx,iPheno),pheno_predict(test_idx,iPheno),1);
            catch
                fold_auc(iPheno,iFold) = NaN;
            end
        end
    end
end

mean_auc = mean(fold_auc,2);
std_auc = std(fold_auc,[],2);
ts_auc = tinv(0.975,nFold-1);
ci_auc = ts_auc.*std_auc./sqrt(nFold);

results.pheno_predict = pheno_predict;
results.pheno_residualized = pheno_residualized;
results.fold_auc = fold_auc;
results.mean_auc = mean_auc;
results.std_auc = std_auc;
results.ts_auc = ts_auc;
results.ci_auc = ci_auc;
results.folds = folds;
