function output_rois = mc_fsaverage5_aseg_rois(value)

structures = {
    8, 'Left-Cerebellum-Cortex'
    10, 'Left-Thalamus'
    11, 'Left-Caudate'
    12, 'Left-Putamen'
    13, 'Left-Pallidum'
    16, 'Brain-Stem'
    17, 'Left-Hippocampus'
    18, 'Left-Amygdala'
    26, 'Left-Accumbens-area'
    28, 'Left-VentralDC'
    47, 'Right-Cerebellum-Cortex'
    49, 'Right-Thalamus'
    50, 'Right-Caudate'
    51, 'Right-Putamen'
    52, 'Right-Pallidum'
    53, 'Right-Hippocampus'
    54, 'Right-Amygdala'
    58, 'Right-Accumbens-area'
    60, 'Right-VentralDC'
};

output_rois = zeros(numel(value),1);
value = reshape(value,numel(value),1);

for i = 1:size(structures,1)
    idx = value==structures{i,1};
    output_rois(idx) = i;
end
