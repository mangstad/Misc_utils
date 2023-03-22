function [cmap,labels] = mc_cmap_redblue(bins)
Networks = {
    'Somatomotor Hand'        '00FFFF';
    'Somatomotor Face'        'FF8000';
    'Cingulo-opercular'       '800080';
    'Auditory'                'FF00FF';
    'Default'                 'FF0000';
    'Unused'                  '808080';
    'Visual'                  '0000FF';
    'Fronto-parietal'         'FFFF00';
    'Salience'                '000000';
    'Subcortical'             'B43C00';
    'Ventral Attention'       '008080';
    'Dorsal Attention'        '00FF00';
    'Cerebellum'              '80FFFF';
    'Uncertain'               'AAAAAA';
    'Cingulo-Parietal'        '804080';
    'Retrosplenial Temporal'  '808080';
    };

rgb = reshape(sscanf(cell2mat(Networks(:,2)).','%2x'),3,[]).'/255;

cmap = rgb
labels = Networks(:,1);