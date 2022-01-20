function [es,pes] = mc_etasq(y,x,groups)
%add intercept if missing, and adjust groups
if (std(x(:,1)) == 0)
    x = [ones(size(x,1),1) x];
    groups = [1 groups];
end

m = fitlm(x,y);

u = unique(groups);
for i = 1:numel(u)
    x0 = x;
    x0(:,groups==u(i)) = [];
    m0 = fitlm(x0,y);
    sst = sum((y-mean(y)).^2);
    %who decided it was a good idea to refer to things as sum of square
    %residuals/error and sum of square effect/regression?
    sse = sum(m.Residuals.Raw.^2);
    ssr = sum((m.Fitted - m0.Fitted).^2);
    es(i) = ssr/sst;
    pes(i) = ssr/(ssr+sse);
end
