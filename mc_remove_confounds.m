function clean_data = mc_remove_confounds(data,confounds,varargin)

remove = 1:(size(confounds,2)+1);

if (nargin>2)
    remove = varargin{1};
end

X = [ones(size(data,1),1) confounds];
b = pinv(X'*X)*X'*data;

clean_data = data - X(:,remove)*b(remove,:);
