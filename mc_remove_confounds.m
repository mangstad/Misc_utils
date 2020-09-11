function clean_data = mc_remove_confounds(data,confounds)

X = [ones(size(confounds,1),1) confounds];
b = pinv(X'*X)*X'*data;

clean_data = data - X*b;
