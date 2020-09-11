function data = mc_load_datafile(File)
%will try to load data based on file type:
%1 = text file matrix (.txt)
%2 = matrices (.mat)
%3 = resting state connectome files from ConnTool (_corr.mat)
%4 = resting state timecourses from ConnTool (_roiTC.mat)
%5 = nifti images (.nii)
%6 = cifti images (.dtseries.nii, .dscalar.nii, etc)

[p,f,e] = fileparts(File);

filetype = 0;
for i = 1:1
    switch e
        case '.txt'
            filetype = 1;
        case '.mat'
            if (~isempty(regexp(f,'_corr$')))
                filetype = 3;
                break;
            end
            if (~isempty(regexp(f,'_roiTC$')))
                filetype = 4;
                break;
            end
            filetype = 2;
        case '.nii'
            if (~isempty(regexp(f,'\.d(tseries)|(label)|(scalar)$')))
                filetype = 6;
                break;
            end
            filetype = 5;
        otherwise
            error('Unsupported file type');
    end
end

switch filetype
    case 1
        data = load(File);
    case 2
        data = load(File);
    case 3
        temp = load(File);
        data = temp.rMatrix;
    case 4
        temp = load(File);
        data = corr(temp.roiTC);
    case 5
        temp = nifti(File);
        data = temp.dat(:,:,:,:);
        s = size(data);
        d = numel(s);
        if (d==4 && s(4)>1)
            data = permute(data,[d 1:(d-1)]);
        end
    case 6
        temp = ciftiopen(File,'/net/parasite/HCP/Scripts/slab/cifti/workbench/bin_linux64/wb_command');
        data = temp.cdata;
        s = size(data);
        d = numel(s);
        if (d==2 && s(2)>1)
            data = permute(data,[d 1:(d-1)]);
        end
end

