function mc_plot_connectome(flatmat,nets,varargin)
    invert = 0;
    if (nargin>2)
        invert = varargin{1};
    end
    h = [];
    if (nargin>3)
        h = varargin{2};
    end
    %simple plotting function to plot a sorted connectome
    [snets,sidx] = sort(nets);
    mat = mc_unflatten_upper_triangle(flatmat);
    mat = mat + mat';
    smat = mat(sidx,sidx);
    if (isempty(h))
        h = figure;
    end
    imagesc(smat);
    %set(gca,'YDir','normal')
    axis square;
    %set colorbar
    cmx = max(abs(smat(:)));
    caxis([-cmx cmx]);
    zo = linspace(0,1,33)';
    oz = flipud(zo);
    zo = zo(1:end-1);
    oz = oz(1:end-1);
    oo = ones(size(zo));
    cmap = [zo zo oo;oo oz oz];
    if (invert)
        cmap = flipud(cmap);
    end
    
    colormap(cmap)
    
    %draw network boundary lines
    hold on;
    u = unique(snets);
    boundaries = cumsum(crosstab(snets));
    boundaries = boundaries(1:end-1);
    limits = axis;
    mn = limits(1);
    mx = limits(2);
    for i = 1:numel(boundaries)
        plot([mn mx],[boundaries(i) boundaries(i)],'k-');
        plot([boundaries(i) boundaries(i)],[mn mx],'k-');
    end
    hold off;
    