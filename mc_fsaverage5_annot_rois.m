function output_rois = mc_fsaverage5_annot_rois(label,ctab)

structures = {
    %'unknown'                 
    'bankssts'                
    'caudalanteriorcingulate' 
    'caudalmiddlefrontal'     
    %'corpuscallosum'          
    'cuneus'                  
    'entorhinal'              
    'fusiform'                
    'inferiorparietal'        
    'inferiortemporal'        
    'isthmuscingulate'        
    'lateraloccipital'        
    'lateralorbitofrontal'    
    'lingual'                 
    'medialorbitofrontal'     
    'middletemporal'          
    'parahippocampal'         
    'paracentral'             
    'parsopercularis'         
    'parsorbitalis'           
    'parstriangularis'        
    'pericalcarine'           
    'postcentral'             
    'posteriorcingulate'      
    'precentral'              
    'precuneus'               
    'rostralanteriorcingulate'
    'rostralmiddlefrontal'    
    'superiorfrontal'         
    'superiorparietal'        
    'superiortemporal'        
    'supramarginal'           
    'frontalpole'             
    'temporalpole'            
    'transversetemporal'      
    'insula'                  
    };

output_rois = zeros(numel(label),1);

for i = 1:numel(structures)
    ctabi = strmatch(structures{i},char(ctab.struct_names));
    code = ctab.table(ctabi,5);
    idx = find(label==code);
    %offset = 0;
    %if (strcmp(hemi,'rh'))
    %    offset = numel(structures);
    %end
    output_rois(idx) = i;%i+offset;
    
end