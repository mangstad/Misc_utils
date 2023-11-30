function [confounds, stats, fd]= fmriprep_getconfounds(filepath,NPC,FDthresh,DetrendOrder,IncludeAroma,IncludeCensor,Trim)
    [path,file,ext] = fileparts(filepath);
    %path = '/nfs/locker/dads-abcd/fmriprep/derivatives/NDARINV0D4C1R8X/sub-NDARINV0D4C1R8X/ses-baselineYear1Arm1/func/';
    %file = 'sub-NDARINV0D4C1R8X_ses-baselineYear1Arm1_task-rest_run-01_desc-confounds_regressors.tsv';
    
    dat = readtable(fullfile(path,[file ext]),'Delimiter','\t','FileType','text','TreatAsEmpty','n/a');
    json = jsondecode(fileread(fullfile(path,[file '.json'])));
    
    f = fields(json);
    %idx = ~cellfun('isempty',strfind(f,'a_comp_cor'));
    idx = contains(f,'a_comp_cor');
    f = f(idx);
    wmv = [];
    wmn = {};
    csfv = [];
    csfn = {};
    for i = 1:numel(f)
        %fprintf(1,'%d\n',i);
        fi = f{i};
        if (strcmp(json.(fi).Method,'aCompCor'))
            switch json.(fi).Mask
                case 'WM'
                    wmv(end+1) = json.(fi).VarianceExplained;
                    wmn{end+1} = fi;
                case 'CSF'
                    csfv(end+1) = json.(fi).VarianceExplained;
                    csfn{end+1} = fi;
            end
        end
    end
    
    [~,i] = sort(wmv,'descend');
    wmns = wmn(i);
    [~,i] = sort(csfv,'descend');
    csfns = csfn(i);
    
    NPCw = min(NPC,numel(wmns)); %adjust number to number present if fewer than NPC
    NPCc = min(NPC,numel(csfns));

    wmn = wmns(1:NPCw);
    csfn = csfns(1:NPCc);
    compcor_regressors = zeros(size(dat,1),NPC*2);
    for i = 1:NPCw
        compcor_regressors(:,i) = dat.(wmn{i});
    end
    for i = 1:NPCc
        compcor_regressors(:,NPCw+i) = dat.(csfn{i});
    end
    
    %AROMA
    f = fields(dat);
    idx = contains(f,'aroma_motion');
    aroma_regressors = table2array(dat(:,idx));
    
    %motion
    f = fields(dat);
    %idx = ~cellfun('isempty',strfind(f,'trans_')) | ~cellfun('isempty',strfind(f,'rot_'));
    idx = contains(f,'trans_') | contains(f,'rot_');
    idxoT = contains(f,'trans_') & ~contains(f,'derivative') & ~contains(f,'power2');
    idxoR = contains(f,'rot_') & ~contains(f,'derivative') & ~contains(f,'power2');
    idxdT = contains(f,'trans_') & contains(f,'derivative') & ~contains(f,'power2');
    idxdR = contains(f,'rot_') & contains(f,'derivative') & ~contains(f,'power2');
    idxqT = contains(f,'trans_') & ~contains(f,'derivative') & contains(f,'power2');
    idxqR = contains(f,'rot_') & ~contains(f,'derivative') & contains(f,'power2');
    idxqdT = contains(f,'trans_') & contains(f,'derivative') & contains(f,'power2');
    idxqdR = contains(f,'rot_') & contains(f,'derivative') & contains(f,'power2');
    
    motion_regressors = table2array([dat(:,idxoT) dat(:,idxoR) dat(:,idxdT) dat(:,idxdR) dat(:,idxqT) dat(:,idxqR) dat(:,idxqdT) dat(:,idxqdR)]);
    
    fd = dat.framewise_displacement;
    censor = fd>FDthresh;
    
    y = censor;
    newy = zeros(numel(y), sum(y));
    newy(sub2ind(size(newy), find(y)',1:sum(y))) = 1;
    censor = newy;
    
    %cosine basis
    f = fields(dat);
    %idx = ~cellfun('isempty',strfind(f,'cosine'));
    idx = contains(f,'cosine');
    cosine_regressors = table2array(dat(:,idx));
    
    %detrend with legendre polynomials as 3dDetrend
    n = [1:size(dat,1)]';
    detrend = [n (1/2)*(3*n.^2-1) (1/2)*(5*n.^3 - 3*n) (1/8)*(35*n.^4 - 30*n.^2 + 3) ...
        (1/8)*(63*n.^5 - 70*n.^3 + 15*n) (1/16)*(231*n.^6 - 315*n.^4 + 105*n.^2 - 5)];
    detrend = detrend(:,1:DetrendOrder);
    
    if (IncludeAroma==0)
        aroma_regressors = [];
    end
    
    if (IncludeCensor==0)
        censor = [];
    end
    
    trimcensor = [];
    if (Trim>0)
        trimcensor = eye(size(dat,1),Trim);
    end
    
    confounds = [ones(size(dat,1),1) detrend motion_regressors compcor_regressors cosine_regressors censor aroma_regressors trimcensor];
    
    %replace NaN with column nanmean
    nm = repmat(nanmean(confounds),size(dat,1),1);
    confounds(isnan(confounds)) = nm(isnan(confounds));
    stats = [sum(censor(:)) size(confounds,1) size(confounds,2) nanmean(dat.framewise_displacement)];
    