function [exp,b,C,v,r2] = VarExpl(Y,X,varargin)

    X = [ones(size(X,1),1) X];

    if (nargin > 2)
        b = varargin{1};
    else
        b = pinv(X)*Y;
    end
    Yhat = X*b;
    
    v = var(Y);
    Y=bsxfun(@minus,Y,mean(Y,1)); %%% zero-mean
    Yhat=bsxfun(@minus,Yhat,mean(Yhat,1)); %%% zero-mean
    Y=bsxfun(@times,Y,1./sqrt(sum(Y.^2,1))); %% L2-normalization
    Yhat=bsxfun(@times,Yhat,1./sqrt(sum(Yhat.^2,1))); %% L2-normalization
    C=sum(Y.*Yhat,1); %% correlation
    r2 = C.^2;
    %r2 = diag(corr(X,Xhat)).^2;
    exp = 100*sum(r2.*(v/sum(v)));
    
    %r = corr(X,Y);
    %r2 = r.^2;
    %v = var(X);
    %vr = repmat(v',1,size(Y,2));
    %exp = sum(sum((vr.*r2)./sum(v)));
    