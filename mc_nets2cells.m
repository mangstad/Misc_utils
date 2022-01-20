function cells = mc_nets2cells(nets,varargin)

dropzero = 0;
if (nargin>1)
    dropzero = varargin{1};
end

u = unique(nets);
if (dropzero)
    u(u==0) = [];
end

cells_square = zeros(numel(nets));
cellidx = 1;
for i = 1:numel(u)
    for j = i:numel(u)
        cells_square(nets==u(i),nets==u(j)) = cellidx;
        cells_square(nets==u(j),nets==u(i)) = cellidx;
        cellidx = cellidx + 1;
    end
end

cellidx = cellidx - 1;

cells = mc_flatten_upper_triangle(cells_square);

