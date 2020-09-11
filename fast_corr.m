function r = fast_corr(a,b)
    a=bsxfun(@minus,a,mean(a,1)); %%% zero-mean
    b=bsxfun(@minus,b,mean(b,1)); %%% zero-mean
    a=bsxfun(@times,a,1./sqrt(sum(a.^2,1))); %% L2-normalization
    b=bsxfun(@times,b,1./sqrt(sum(b.^2,1))); %% L2-normalization
    r=sum(a.*b,1); %% correlation