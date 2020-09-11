function a = plot_jica_component(i,flip,ThresholdType,Threshold,nets,varargin)
    cdzfx = flip;
    
    cCon = i;
    cConSign = sign(cCon);
    
    title = '';
    if (nargin>5)
        title = varargin{1};
    end
    
    subset = 1:12;
    if (nargin>6)
        subset = varargin{2};
    end
    
    h = [];
    if (nargin>7)
        h = varargin{3};
    end
    
    switch (ThresholdType)
        case 0
            cConThresh = cCon;
            cConThresh(abs(cCon) < Threshold) = 0;
        case 1
            cConZ = zscore(cCon);
            cConThresh = cConZ;
            cConThresh(abs(cConZ) < Threshold) = 0;
        case 2
            pct = prctile(abs(cCon),100-Threshold);
            cConThresh = cCon;
            cConThresh(abs(cConThresh)<pct) = 0;
    end
    
    cConThresh = cConThresh * cdzfx;
    
    clear a   
    if(~isempty(h))
        a.h = h;
    end
    a.mediator.NetSubset =  subset;
    a.tvalues = zeros(size(cConThresh));
    a.tvalues(sign(cConThresh)== 0) = 1;
    a.tvalues(sign(cConThresh)==+1) = 2;
    a.tvalues(sign(cConThresh)==-1) = 3;
    a.bvalues = a.tvalues; % just satisfy the function
    a.bvalues = cConThresh;
    a.NetworkLabels = nets;
    a.dotenable = 0;
    a = mc_Network_FeatRestruct(a);
    %a.title = ['Component TakGraph ' title];
    a.title = title;
    a = mc_TakGraph_plot(a);
end
