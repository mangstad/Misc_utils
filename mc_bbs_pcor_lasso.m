function results = mc_bbs(featuremat,pheno,nuisance,folds,NumComp,varargin)

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
addParameter(p,'NoComponents',0);
addParameter(p,'Logistic',0); %currently not working
addParameter(p,'Seed',-1);

parse(p,varargin{:});

n = size(featuremat,1);

nFold = numel(unique(folds));

NumComps = numel(NumComp);

testpheno = p.Results.TestPheno;

pcadone = 0;
if (~isempty(p.Results.Scores))
    Aa = p.Results.Scores;
    pcadone = 1;
    results.Aa = Aa;
end
results.components = [];
if (~isempty(p.Results.Components))
    components = p.Results.Components;
    results.components = components;
end
LOSOPheno = p.Results.LOSOPheno;

NoComponents = p.Results.NoComponents;

Logistic = p.Results.Logistic;

Seed = p.Results.Seed;

clear p;

good = ~any(isnan(pheno),2) & ~any(isnan(nuisance),2);

nFold = numel(unique(folds));
if (pcadone==0) 
clear Aa components;
    for iFold = 1:nFold
        tic
        fprintf(1,'.');
        %find train and test data for this fold
        test_idx = folds==iFold;% & good;
        train_idx = ~test_idx;% & good;
        
        %reduce the training data
        coeff = pca(featuremat(train_idx,:));
        components{iFold} = coeff';

        %mean center train, and mean center test with train means
        mu = mean(featuremat(train_idx,:));
        x = bsxfun(@minus,featuremat,mu);

        %calculate expressions for each subject for train and test
        Aa{iFold} = (pinv(components{iFold}')*x')';
        toc
    end
    results.Aa = Aa;
    if (NoComponents==0)
        results.components = components;
    end
end

if (Logistic==1)
    %results = mc_bbs_logistic(featuremat,pheno,nuisance,folds,NumComp,'Scores',results.Aa,'Components',results.components,'LOSOPheno',LOSOPheno,'TestPheno',testpheno);
    error('Logistic currently not supported with partial corr method.');
    return;
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

    %if NumComps>1 then do nested 5-fold
    if (NumComps>1)
        if (Seed==-1)
            s = rng('shuffle');
            Seed = s.Seed;
        end
        rng(Seed);
        results.Seed = Seed;
        nestfold = randsample(5,n_trainf,1);
        %nestAa = [];
        %for iNest = 1:10
        %    nestAa{iNest} = Aa{iFold}(train_idx,:);
        %end
        nestresults = zeros(size(pheno,2),NumComps);
        firstrun = mc_bbs_pcor(featuremat(train_idx,:),pheno(train_idx,:),nuisance(train_idx,:),nestfold,1,'LOSOPheno',LOSOPheno);
        for iNest = 1:NumComps
            tempresults = mc_bbs_pcor(featuremat(train_idx,:),pheno(train_idx,:),nuisance(train_idx,:),nestfold,NumComp(iNest),'Scores',firstrun.Aa,'LOSOPheno',LOSOPheno);
            nestresults(:,iNest) = tempresults.mean_corr;
        end
        [bestresults,bestcomps] = max(nestresults,[],2);
        %calculate 1se of results over fold
        onese = std(tempresults.fold_corr)/sqrt(numel(tempresults.fold_corr));
        min1se = bestresults-onese;
        oneseidx = find(nestresults>min1se);
        onesecomps = NumComp(oneseidx(1))';
        bestcomps = NumComp(bestcomps)';
        results.bestcomps(:,iFold) = bestcomps;
        results.onesecomps(:,iFold) = onesecomps;
        
        Abig = Aa{iFold}(train_idx,1:max(NumComp));
        Abig_test = Aa{iFold}(test_idx,1:max(NumComp));
        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            kk = bestcomps(iPheno);
            
            %partial corr residualize in train
            X = [ones(n_trainf,1) nuisance(train_idx,:)];
            b1 = pinv(X'*X)*X'*pheno(train_idx,iPheno);
            b2 = pinv(X'*X)*X'*Abig(:,1:kk);
            ytr = pheno(train_idx,iPheno) - X*b1;
            xtr = Abig(:,1:kk) - X*b2;
            
            Xte = [ones(n_testf,1) nuisance(test_idx,:)];
            yte = pheno(test_idx,iPheno) - Xte*b1;
            xte = Abig_test(:,1:kk) - Xte*b2;
            
            X = [ones(n_trainf,1) xtr];
            b = pinv(X'*X)*X'*ytr;
            
            %X = [ones(n_trainf,1) Abig(:,1:kk) nuisance(train_idx,:)];
            %b = pinv(X'*X)*X'*pheno(train_idx,iPheno);

            %pheno_predict(test_idx,iPheno) = Abig_test(:,1:kk)*b(2:(kk+1));
            pheno_predict(test_idx,iPheno) = xte(:,1:kk)*b(2:(kk+1));
            
            if (isempty(testpheno))
                %pheno_residualized(test_idx,iPheno) = pheno(test_idx,iPheno) - [ones(n_testf,1) nuisance(test_idx,:)]*b([1 (kk+2):end]);
                pheno_residualized(test_idx,iPheno) = yte;
            else
                pheno_residualized(test_idx,iPheno) = testpheno(test_idx,iPheno) - [ones(n_testf,1) nuisance(test_idx,:)]*b1;
            end

        end
    else
        Abig = Aa{iFold}(train_idx,1:k);
        Abig_test = Aa{iFold}(test_idx,1:k);
        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            
            X = [ones(n_trainf,1) nuisance(train_idx,:)];
            b1 = pinv(X'*X)*X'*pheno(train_idx,iPheno);
            b2 = pinv(X'*X)*X'*Abig;
            Abig_r = Abig - X*b2;
            pheno_r = pheno(train_idx,iPheno) - X*b1;
            
            Xt = [ones(n_testf,1) nuisance(test_idx,:)];
            Abig_test_r = Abig_test - Xt*b2;
            pheno_test_r = pheno(test_idx,iPheno) - Xt*b1;
            
            Xm = [ones(n_trainf,1) Abig_r];
            bm = pinv(Xm'*Xm)*Xm'*pheno_r;
            
            Xmt = [ones(n_testf,1) Abig_test_r];
            
            pheno_predict(test_idx,iPheno) = Xmt*bm;
            
            if (isempty(testpheno))
                pheno_residualized(test_idx,iPheno) = pheno_test_r;
            else
                pheno_residualized(test_idx,iPheno) = testpheno(test_idx,iPheno) - Xt*b1;
            end
        end
    end
end

%check correlation between actual and predicted phenotypes
fold_corr = zeros(size(pheno,2),nFold);
fold_mse = zeros(size(pheno,2),nFold);
fold_nmse = zeros(size(pheno,2),nFold);
fold_r2cv = zeros(size(pheno,2),nFold);
fold_petasq = zeros(size(pheno,2),nFold);

if (LOSOPheno)
    fold_corr = zeros(1,nFold);
    fold_mse = zeros(1,nFold);
    fold_nmse = zeros(1,nFold);
    fold_r2cv = zeros(1,nFold);
    fold_petasq = zeros(1,nFold);
end

for iFold = 1:nFold
    test_idx = folds==iFold & good;
    train_idx = ~test_idx & good;
    if (LOSOPheno)
        fold_corr(:,iFold) = diag(corr(pheno_predict(test_idx,iFold),pheno_residualized(test_idx,iFold)));
        fold_mse(:,iFold) = mean((pheno_residualized(test_idx,iFold)-pheno_predict(test_idx,iFold)).^2);
        fold_nmse(:,iFold) = fold_mse(:,iFold)./mean(bsxfun(@minus,pheno_residualized(test_idx,iFold),mean(pheno_residualized(train_idx,iFold))).^2)';
        fold_r2cv(:,iFold) = 1-fold_nmse(:,iFold);
        [~,pe] = mc_etasq(pheno_residualized(test_idx,iFold),pheno_predict(test_idx,iFold),1);
        fold_petasq(:,iFold) = pe; 
    else
        fold_corr(:,iFold) = diag(corr(pheno_predict(test_idx,:),pheno_residualized(test_idx,:)));
        fold_mse(:,iFold) = mean((pheno_residualized(test_idx,:)-pheno_predict(test_idx,:)).^2);
        fold_nmse(:,iFold) = fold_mse(:,iFold)./mean(bsxfun(@minus,pheno_residualized(test_idx,:),mean(pheno_residualized(train_idx,:))).^2)';
        fold_r2cv(:,iFold) = 1-fold_nmse(:,iFold);
        for iPheno = 1:size(pheno,2)
            [~,pe] = mc_etasq(pheno_residualized(test_idx,iPheno),pheno_predict(test_idx,iPheno),1);
            fold_petasq(iPheno,iFold) = pe;
        end
        
    end
end

mean_corr = mc_FisherZ(mean(mc_FisherZ(fold_corr),2),1);
std_corr = mc_FisherZ(std(mc_FisherZ(fold_corr),[],2),1);
ts_corr = tinv(0.975,nFold-1);
ci_corr = ts_corr.*std_corr./sqrt(nFold);

results.fold_mse = fold_mse;
results.fold_nmse = fold_nmse;
results.fold_r2cv = fold_r2cv;
results.fold_petasq = fold_petasq;

results.mean_mse = mean(fold_mse,2);
results.mean_nmse = mean(fold_nmse,2);
results.mean_r2cv = mean(fold_r2cv,2);
results.mean_petasq = mean(fold_petasq,2);

results.pheno_predict = pheno_predict;
results.pheno_residualized = pheno_residualized;
results.fold_corr = fold_corr;
results.mean_corr = mean_corr;
results.std_corr = std_corr;
results.ts_corr = ts_corr;
results.ci_corr = ci_corr;
results.folds = folds;
