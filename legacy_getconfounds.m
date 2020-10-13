function [confounds, stats]= legacy_getconfounds(data,filepathM,filepathW,filepathC,TR,HPF,NPC,FDthresh,DetrendOrder,IncludeCensor,Trim)
    motion = load(filepathM);
    wmmask = mc_load_datafile(filepathW);
    csfmask = mc_load_datafile(filepathC);

    wmmask = reshape(wmmask,1,numel(wmmask));
    csfmask = reshape(csfmask,1,numel(csfmask));
    
    [~,wmns] = pca(data(:,wmmask));
    [~,csfns] = pca(data(:,csfmask));
    
    wmn = wmns(:,1:NPC);
    csfn = csfns(:,1:NPC);
    compcor_regressors = [wmn csfn];

    %motion
    dmotion = [repmat(0,1,size(motion,2));diff(motion)];
    
    motion_regressors = [motion dmotion motion.^2 dmotion.^2];
    
    %calculate framewise displacement
    %
    
    tmp = dmotion;
    tmp(:,4:6) = tmp(:,4:6)*50;
    fd = sum(abs(tmp),2);
    censor = fd>FDthresh;
    
    y = censor;
    newy = zeros(numel(y), sum(y));
    newy(sub2ind(size(newy), find(y)',1:sum(y))) = 1;
    censor = newy;
    
    %cosine basis
    frametimes = linspace(0,TR*size(data,1)-TR,size(data,1));
    cosine_regressors = cosine_basis(1/(2*TR),HPF,frametimes);
    
    %detrend with legendre polynomials as 3dDetrend
    n = [1:size(dat,1)]';
    detrend = [n (1/2)*(3*n.^2-1) (1/2)*(5*n.^3 - 3*n) (1/8)*(35*n.^4 - 30*n.^2 + 3) ...
        (1/8)*(63*n.^5 - 70*n.^3 + 15*n) (1/16)*(231*n.^6 - 315*n.^4 + 105*n.^2 - 5)];
    detrend = detrend(:,1:DetrendOrder);
    
    if (IncludeCensor==0)
        censor = [];
    end
    
    trimcensor = [];
    if (Trim>0)
        trimcensor = eye(size(dat,1),Trim);
    end
    
    confounds = [ones(size(dat,1),1) detrend motion_regressors compcor_regressors cosine_regressors censor trimcensor];
    
    %replace NaN with column nanmean
    nm = repmat(nanmean(confounds),size(dat,1),1);
    confounds(isnan(confounds)) = nm(isnan(confounds));
    stats = [sum(censor(:)) size(confounds,1) size(confounds,2) nanmean(fd)];
    