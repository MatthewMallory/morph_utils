import os
from tqdm import tqdm
import argschema as ags
import numpy as np
from morph_utils.ccf import open_ccf_annotation
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from allensdk.core.reference_space_cache import ReferenceSpaceCache
try:   
    from skeleton_keys import full_morph
    
except ImportError:
    msg = """
    Required module (skeleton_keys.full_morph) is not installed. It's possible you have skeleton_keys installed
    but not the correct branch/version. As of 12/22/23 the full_morph branch has not been merged into the main 
    branch of skeleton_keys so check the full_morph-MM-edits branch for the full_morph features used in this code
    """
    print(msg)
    
    
class IO_Schema(ags.ArgSchema):
    input_ccf_swc_directory = ags.fields.InputDir(description='directory with micron resolution ccf registered files')
    output_ccf_swc_directory = ags.fields.OutputDir(description="output directory for swc files")
    closest_surface_voxel_file = ags.fields.InputFile(description="path to closest_surface_voxel_file")
    surface_paths_file = ags.fields.InputFile(description="path to surface_paths_file")

def main(input_ccf_swc_directory, 
         output_ccf_swc_directory, 
         closest_surface_voxel_file, 
         surface_paths_file,
         **kwargs):
    
    atlas_volume = open_ccf_annotation(with_nrrd=True)
    
    reference_space_key = 'annotation/ccf_2017'
    resolution = 10
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    tree = rspc.get_structure_tree(structure_graph_id=1)
    acronym_map = rspc.get_reference_space().structure_tree.get_id_acronym_map()
    
    isocortex_struct_id = acronym_map['Isocortex']
    for swc_fn in tqdm([f for f in os.listdir(input_ccf_swc_directory) if ".swc" in f]):

        swc_pth = os.path.join(input_ccf_swc_directory, swc_fn)
        output_swc_pth = os.path.join(output_ccf_swc_directory, swc_fn)
        
        morph = morphology_from_swc(swc_pth)
        morph_soma = morph.get_soma()
        soma_coords = np.array([morph_soma['x'], morph_soma['y'], morph_soma['z']])

        out_of_cortex_bool, nearest_cortex_coord = full_morph.check_coord_out_of_cortex(soma_coords,
                                                                                isocortex_struct_id,
                                                                                atlas_volume,
                                                                                closest_surface_voxel_file,
                                                                                surface_paths_file,
                                                                                tree)

        if out_of_cortex_bool:
            soma_coords = nearest_cortex_coord
        cropped_morph = full_morph.local_crop_cortical_morphology(morph, 
                                                        soma_coords, 
                                                        closest_surface_voxel_file, 
                                                        surface_paths_file,
                                                        threshold=500)
        
        morphology_to_swc(cropped_morph, output_swc_pth)

def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)

if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
