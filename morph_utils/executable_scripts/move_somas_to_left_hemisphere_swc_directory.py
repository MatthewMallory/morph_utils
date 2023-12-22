import os
from tqdm import tqdm
import pandas as pd
import argschema as ags
from morph_utils.ccf import move_soma_to_left_hemisphere
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc


class IO_Schema(ags.ArgSchema):
    input_ccf_swc_directory = ags.fields.InputDir(description='directory with micron resolution ccf registered files')
    output_ccf_swc_directory = ags.fields.OutputFile(description="output directory for swc files")
    resolution = ags.fields.Int(default=10, description="Optional. ccf resolution (micron/pixel")
    volume_shape = ags.fields.List(ags.fields.Int, default=[1320, 800, 1140], description = "Optional. Size of input annotation")

def main(input_ccf_swc_directory, 
         output_ccf_swc_directory, 
         resolution, 
         volume_shape,
         **kwargs):
    
    z_size = resolution * volume_shape[2]
    z_midline = z_size / 2

    for swc_fn in tqdm([f for f in os.listdir(input_ccf_swc_directory) if ".swc" in f]):

        swc_pth = os.path.join(input_ccf_swc_directory, swc_fn)
        output_swc_pth = os.path.join(output_ccf_swc_directory, swc_fn)
        
        morph = morphology_from_swc(swc_pth)
        morph = move_soma_to_left_hemisphere(morph, resolution, volume_shape, z_midline)
        morphology_to_swc(morph, output_swc_pth)

def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)

if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
