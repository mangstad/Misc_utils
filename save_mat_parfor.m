function save_mat_parfor(path,mat,varargin)
    ascii = 1;
    if (nargin>2)
        ascii = varargin{1};
    end
    
    if (ascii==1)
        save(path,'mat','-ASCII');
    else
        save(path,'mat','-v7.3');
    end
    