function cmap = mc_cmap_redblue(bins)

zo = linspace(0,1,bins)';
oz = flipud(zo);
o = ones(bins,1);

cmap = [[zo zo o];[1 1 1];[o oz oz]];
