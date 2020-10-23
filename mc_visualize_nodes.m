function image = mc_visualize_nodes(values,atlas,output)

%need code to deal with cifti vs nifti, but for now just cifti

addpath /net/parasite/HCP/Scripts/slab/cifti
addpath /net/parasite/HCP/Scripts/slab/cifti/gifti-1.6

nn = ciftiopen(atlas,'/net/parasite/HCP/Scripts/slab/cifti/workbench/bin_linux64/wb_command');
nnout = nn;
nodes = nn.cdata;
data = 0*nodes;

for i = 1:numel(values)
    data(nodes==i) = values(i);
end

nnout.cdata = data;
ciftisavereset(nnout,output,'/net/parasite/HCP/Scripts/slab/cifti/workbench/bin_linux64/wb_command');

image = nnout;
