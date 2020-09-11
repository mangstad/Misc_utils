function [confounds, stats]= fmriprep_getconfounds(filepath,NPC,FDthresh,DoCensor,Trim)
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
    
    wmn = wmns(1:NPC);
    csfn = csfns(1:NPC);
    compcor_regressors = zeros(size(dat,1),NPC*2);
    for i = 1:NPC
        compcor_regressors(:,i) = dat.(wmn{i});
        compcor_regressors(:,NPC+i) = dat.(csfn{i});
    end
    
    %motion
    f = fields(dat);
    %idx = ~cellfun('isempty',strfind(f,'trans_')) | ~cellfun('isempty',strfind(f,'rot_'));
    idx = contains(f,'trans_') | contains(f,'rot_');
    motion_regressors = table2array(dat(:,idx));
    
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
    
    trimcensor = eye(size(dat,1),Trim);
    
    %confounds = [ones(size(dat,1),1) compcor_regressors motion_regressors cosine_regressors trimcensor];
    confounds = [compcor_regressors motion_regressors cosine_regressors trimcensor];
    if (DoCensor)
        confounds = [confounds censor];
    end

    %replace NaN with column nanmean
    nm = repmat(nanmean(confounds),size(dat,1),1);
    confounds(isnan(confounds)) = nm(isnan(confounds));
    stats = [sum(censor(:)) size(confounds,1) size(confounds,2) nanmean(dat.framewise_displacement)];
    