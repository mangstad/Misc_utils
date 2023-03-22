function a = mc_nca_plot(results,thresh,col,nets,cm,fdr)

if (size(nets,1)>1)
    nets = nets';
end
%a = plot_jica_component(results.stats.t(col,:),1,0,0,nets,'',[1:max(nets)]);
% 
% a.tvalues = results.stats.t(col,:).*(abs(results.stats.t(col,:))>results.tthresh(thresh));
% a.bvalues = results.stats.b(col,:).*(abs(results.stats.t(col,:))>results.tthresh(thresh));
% a.values = a.tvalues;
% a.NetworkLabels = nets;
% a.dotenable = 1;
% a.title = '';
% 
% a.DotDilateMat = [1 0; -1 0; 0 1; 0 -1; % cross
%                    -1 1; 1 1; -1 -1; 1 -1]; %fill out square
%                    % -2 0; 0 2; 2 0; 0 -2]; % cross around square
% a.colormap = [1 1 1; % make 1 white
%               1 0 0; % make 2 red
%               0 0 1; % make 3 blue
%               1 1 0; % make 4 yellow (blended)
%               1 0 0; % make 5 dark red
%               0 0 1; % make 6 dark blue
%                     ];
% a.mediator.NetSubset = 1:max(nets);
%a.mediator.NetSubset =  subset;
tvals = results.stats.t(col,:);
tvals(isnan(tvals)) = 0;
tvals(isinf(tvals)) = 0;
tvals(abs(tvals) < results.tthresh(thresh)) = 0;
a.tvalues = zeros(size(tvals));
a.tvalues(sign(tvals)== 0) = 1;
a.tvalues(sign(tvals)==+1) = 2;
a.tvalues(sign(tvals)==-1) = 3;
a.bvalues = results.stats.t(col,:);
a.NetworkLabels = nets;
a.dotenable = 0;
a.title = 'plot';

a.stats.FDR.NetIn = [1:max(nets)];
a.stats.FDR.Enable = fdr;
a.stats.FDR.rate  = results.FDRthresh;
a.stats.FDR.mode  = 'pdep';
a.stats.FDR.CalcP    = 1;                

a.perms = results.shufcell_count(:,:,:,thresh,col);
a.meanbperms = results.shufcell_mean(:,:,:,thresh,col);
                    
a = mc_Network_FeatRestruct(a);
a = mc_Network_CellCount(a);
%a = mc_TakGraph_plot(a);
a = mc_Network_CellLevelStats(a);
a = mc_TakGraph_CalcShadeColor(a,cm);

if (fdr==0)
    a.shading.shademask = a.stats.rawp<results.FDRthresh(thresh);
end
a = mc_TakGraph_plot(a);
a = mc_TakGraph_AddShading(a,cm);
