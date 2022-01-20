function [data,errors] = mc_load_data(Template,subs,varargin)

flatten = 1;
if (nargin>2)
    flatten = varargin{1};
end

N = length(subs);

Subject = subs{1};
tmp = mc_load_datafile(mc_GenPath(Template));
s = size(tmp);
s2 = [N s];
data = zeros(s2);

if (flatten)
    tmp = reshape(tmp,1,numel(tmp));
    data = reshape(data,N,prod(s));
    data(1,:) = tmp;
else
    %need to generalize this to N-d data
    data(1,:,:,:) = tmp;
end

errors = cell(N,1);

for i = 2:N
    Subject = subs{i};
    fprintf(1,'%d of %d\n',i,N);
    [tmppath,errors{i}] = mc_GenPath(Template);
    if (~isempty(errors{i}))
        continue;
    end
    tmp = mc_load_datafile(tmppath);
    if (flatten)
        tmp = reshape(tmp,1,numel(tmp));
        data(i,:) = tmp;
    else
        data(i,:,:,:) = tmp;
    end
end
