function mtx = mc_threedmap(template,VnewPath,data,roiMNI,varargin)
% MC_THREEDMAP
% Reconstruct 3d nii image with node-wise data.
%
% INPUT
% mask/template - Mask for building your map. This also serves as the template to define voxels sizes/etc.
% VnewPath      - The new nii file path
% data          - nVoxel x 1 (or 1 x nVoxel) vector, the data that being written in
%                 the new nii file.
% roiMNI        - mni coordinates
%
% OPTIONAL INPUT
% connectivity - {1} = single voxel, {7} = faces, {19} = faces + edges, {27} =
% faces + edges + corners (cube), any value not in a cell = radius in
% voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin~=5
    expand = 1;
elseif (iscell(varargin{1}))
    expand = varargin{1}{1};
else
    expand = 0;
end

% Select template file
Vmask = spm_vol(template);
Vnew = Vmask;
Vnew=rmfield(Vnew,'descrip');
Vnew=rmfield(Vnew,'pinfo');

% read mask file
mask  = spm_read_vols(Vmask);

% convert coordinates
roiVoxels = inv(Vnew.mat)*[roiMNI ones(size(roiMNI,1),1)]';
roiVoxels = round(roiVoxels(1:3,:))';

% Initiate output data
Vnew.fname = VnewPath; 
mtx        = zeros(Vnew.dim);

switch expand
    case 0
        radius = varargin{1};
    case 1
        radius = .9;
    case 7
        radius = 1;
    case 19
        radius = 1.5;
    case 27
        radius = 2;
end
  
%[x,y,z] = ndgrid(-1:1);
sz = max(1,ceil(radius));
[x,y,z] = ndgrid(-sz:sz);

kernel = sqrt(x.^2 + y.^2 + z.^2) <=radius; %determine this based on connectivity
modifiers = [x(kernel) y(kernel) z(kernel)];

xlim = Vnew.dim(1);
ylim = Vnew.dim(2);
zlim = Vnew.dim(3);

for i = 1:length(data)
    lidx = [];
    idxx = roiVoxels(i,1);
    idxy = roiVoxels(i,2);
    idxz = roiVoxels(i,3);    
    xyz = repmat([idxx idxy idxz],size(modifiers,1),1);
    xyz = xyz + modifiers;
    lidx = sub2ind(Vnew.dim,xyz(:,1),xyz(:,2),xyz(:,3));
    mtx(lidx) = data(i);
end


% write out results
mtx       = mtx.*mask;
spm_write_vol(Vnew,mtx);

end
