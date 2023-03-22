function results = mc_nca(featuremat,pheno,edgethresh,FDRthresh,nets,varargin)
p = inputParser;
%addRequired(p,'featuremat');
%addRequired(p,'pheno');
%addRequired(p,'edgethresh');
%addRequired(p,'FDRthresh');
%addRequired(p,'nets');
addParameter(p,'nPerm',0);
addParameter(p,'nuisance',[]);
addParameter(p,'verbose',0);
addParameter(p,'Perms',[]);
addParameter(p,'cols',[]);
addParameter(p,'abs',0);

parse(p,varargin{:});
nuisance = p.Results.nuisance;
nPerm = p.Results.nPerm;
Perms = p.Results.Perms;
verbose = p.Results.verbose;
absmean = p.Results.abs;

n = size(featuremat,1);

if (isempty(nuisance))
    nuisance = zeros(n,1);
end

nNets = numel(unique(nets));
[nets2,netsmapping] = mc_correct_folds(nets);

cols = p.Results.cols;
if (isempty(cols))
    cols = 1 + [1:size(pheno,2)];
end

percell_count = zeros(nNets,nNets,numel(edgethresh),numel(cols));
percell_mean = zeros(nNets,nNets,numel(edgethresh),numel(cols));

%run regression to get edgewise p values for effect of interest
X = [ones(size(featuremat,1),1) pheno nuisance];
good = ~any(isnan(X),2);
%tthresh = tinv(1-(edgethresh/2),sum(good)-size(X,2));
statsp = mc_CovariateCorrectionFaster(featuremat(good,:),X(good,:),3,cols,1);

%now calculate cell counts for things above threshold
for iThresh = 1:numel(edgethresh)
    tthresh = tinv(1-(edgethresh(iThresh)/2),sum(good)-size(X,2));
    results.tthresh(iThresh) = tthresh;
    for iP = 1:numel(cols)
        for i = 1:nNets
            for j = i:nNets
                mask = zeros(numel(nets));
                mask(nets2==i,nets2==j) = 1;
                mask(nets2==j,nets2==i) = 1;
                flatmask = mc_flatten_upper_triangle(mask)==1;
                if (sum(flatmask)>0)
                    percell_count(i,j,iThresh,iP) = sum(statsp.t(iP,flatmask)>tthresh);
                    if (absmean)
                        percell_mean(i,j,iThresh,iP) = mean(abs(statsp.t(iP,flatmask)));
                    else
                        percell_mean(i,j,iThresh,iP) = mean(statsp.t(iP,flatmask));
                    end
                end
            end
        end
    end
end



shuf_idx = zeros(n,nPerm);
if (isempty(Perms)) %just generate within-fold permutations
    for iPerm = 1:nPerm
        shuf_idx(:,iPerm) = randperm(sum(good));
        shuf_sign(:,iPerm) = ones(sum(good),1);
    end
else 
    %use provied permutation matrix
    shuf_idx = abs(Perms);
    shuf_sign = sign(Perms);
end

%now get residuals of edges with respect to nuisance
sfm = featuremat(good,:);
sX = X(good,:);
nX = [ones(size(featuremat,1),1) nuisance];
nX = nX(good,:);
b = pinv(nX'*nX)*nX'*sfm;
res = sfm - nX(:,2:end)*b(2:end,:);

%now do perms, freedman and lane style
%that is, first you remove effect of confounds to get residuals
%then you shuffle those residuals
%then you add back true effect of nuisance
%then recalculate your original regression
%permout = zeros(size(featuremat,2),size(pheno,2),nPerm);
nXb = single(nX*b);
res = single(res);
shufcell_count = zeros(nNets,nNets,nPerm,numel(edgethresh),numel(cols));
shufcell_mean = zeros(nNets,nNets,nPerm,numel(edgethresh),numel(cols));
    
for iPerm = 1:nPerm
    if(verbose==1)
        fprintf(1,'%d\n',iPerm);
    end
    shuffle = shuf_idx(:,iPerm);
    sres = res(shuffle,:);
    %%b = a(bsxfun(@plus, ind, 0:size(a,1):numel(a)-1)); %// convert to linear index
    %%sres = res(bsxfun(@plus,shuffle,0:size(res,1):numel(res)-1));
    sres = bsxfun(@times,sres,shuf_sign(:,iPerm)) + nXb;
    shufstats = mc_CovariateCorrectionFaster(sres,sX,3,cols,1);
    %permout(:,:,iPerm) = shufstats.t';
    permout = shufstats.t';

    for iThresh = 1:numel(edgethresh)
    %now calculate cell counts for permuted data
        for iP = 1:numel(cols);
            for i = 1:nNets
                tmp = zeros(nNets,1);
                tmp2 = zeros(nNets,1);
                for j = i:nNets
                    mask = zeros(numel(nets));
                    mask(nets2==i,nets2==j) = 1;
                    mask(nets2==j,nets2==i) = 1;
                    flatmask = mc_flatten_upper_triangle(mask)==1;
                    if (sum(flatmask)>0)
                        tmp(j) = sum(abs(permout(flatmask,iP))>tthresh);
                        if (absmean)
                            tmp2(j) = mean(abs(permout(flatmask,iP)));
                        else
                            tmp2(j) = mean(permout(flatmask,iP));
                        end
                    end
                end
                shufcell_count(i,:,iPerm,iThresh,iP) = tmp;
                shufcell_mean(i,:,iPerm,iThresh,iP) = tmp2;
            end
        end
    end
    
end

%now we'll borrow a bunch of the old NCA code to plot and calculate FDR
%corrected p values for cells
cm = 1; % 1 for count, 2 for mean

results.edgethresh = edgethresh;
results.FDRthresh = FDRthresh;
results.pheno = pheno;
results.stats = statsp;
results.shufcell_count = shufcell_count;
results.shufcell_mean = shufcell_mean;
%results.perms = permout;
results.count = percell_count;
results.mean = percell_mean;

results.netsmapping = netsmapping;


% for pidx = 1:size(pheno,2)
%     %rmpath('/home/slab/users/mangstad/repos/MethodsCore/matlabScripts/Takgraph');
%     tmp = (statsp.t(pidx,:).*(abs(statsp.t(pidx,:))>tthresh));
%     a = plot_jica_component(tmp,1,0,0,nets2','',[1:nNets]);
%     %addpath('/home/slab/users/mangstad/repos/MethodsCore/matlabScripts/Takgraph');
%     a.values = a.tvalues;
%     a = mc_Network_CellCount(mc_Network_FeatRestruct(a));
%     
%     a.perms = shufcell_count(:,:,:,pidx);
%     a.meanbperms = shufcell_mean(:,:,:,pidx);
%     
%     a.stats.FDR.NetIn = [1:nNets];
%     a.stats.FDR.Enable = 1;
%     a.stats.FDR.rate  = FDRthresh;
%     a.stats.FDR.mode  = 'pdep';
%     a.stats.FDR.CalcP    = 1;
%     
%     a = mc_Network_CellLevelStats(a);
%     a = mc_TakGraph_CalcShadeColor(a,cm);
%     a = mc_TakGraph_AddShading(a,cm);
%     results.a(pidx) = a;
% end




