function dummy = mc_nuisance_dummy(nuisance,varargin)

    drop_biggest = 1;
    if (nargin>1)
        drop_biggest = varargin{1};
    end
    
    if (iscellstr(nuisance))
        nans = strcmp(nuisance,'NaN');
    else
        nans = isnan(nuisance); %fix to deal with nans in vector
    end
    
    dummy = mc_Vector2Mask(nuisance);
    dummy(nans,:) = NaN;
    i = nansum(dummy)==0;
    dummy(:,i) = [];
    [~,i] = max(nansum(dummy));
    if (drop_biggest)
        dummy(:,i) = [];
    end
    