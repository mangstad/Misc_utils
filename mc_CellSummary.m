function [cellcount,cellsum,cellsize] = mc_CellSummary(map,nets)
u = unique(nets);

cellsize = [];
cellcount = [];
cellsum = [];

idx = 1;
for i = 1:numel(u)
    for j = i:numel(u)
        ii = u(i);
        jj = u(j);
        mask = zeros(418,418);
        mask(nets==ii,nets==jj) = 1;
        mask(nets==jj,nets==ii) = 1;
        mask = mc_flatten_upper_triangle(mask)==1;
        if (sum(mask)>0)
            cellsize(idx) = sum(mask);
            cellcount(idx) = sum(map(mask)~=0);
            cellsum(idx) = sum(map(mask));
            idx = idx + 1;
        end
    end
end

cellsize = cellsize';
cellcount = cellcount';
cellsum = cellsum';
