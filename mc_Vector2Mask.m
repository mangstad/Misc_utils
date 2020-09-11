function mask = mc_Vector2Mask(nets,varargin)
    scale = 0;
    if (nargin>1)
        scale = varargin{1};
    end
    
    [u,~,uidx] = unique(nets);
    newidx = 1:numel(u);
    newnets = newidx(uidx);
    mask = zeros(numel(nets),numel(newidx));
    mask(sub2ind(size(mask),1:numel(nets),newnets))=1;
    if (scale)
        mask = bsxfun(@rdivide,mask,sum(mask));
    end
    
    
