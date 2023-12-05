import os
import ast
import json
import nrrd
from importlib.resources import files
import pandas as pd
import numpy as np
import SimpleITK as sitk
from neuron_morphology.swc_io import morphology_from_swc
from neuron_morphology.transforms.affine_transform import AffineTransform as aff
from morph_utils.query import get_id_by_name, get_structures, query_pinning_info_cell_locator
from morph_utils.measurements import get_node_spacing

NAME_MAP_FILE = files('morph_utils') / 'data/ccf_structure_name_map.json'
with open(NAME_MAP_FILE, "r") as fn: 
    NAME_MAP = json.load(fn)
NAME_MAP = {int(k):v for k,v in NAME_MAP.items()}

def open_ccf_annotation(with_nrrd, annotation_path=None):
    """
    Open up CCF annotation volume. Use nrrd to open file to get 3-d array, or set with_nrrd to false 
    to open with Sitk. These result in different data structures.

    Args:
        with_nrrd (bool): True if you want to use nrrd to open file, False if you want to use sitk.ReadImage
        annotation_path (str, optional): path to annotation.nrrd file. Defaults to None.

    Returns:
        array: 3d atlas array
    """
    if annotation_path is None: 
        annotation_path =  files('morph_utils') / 'data/annotation_10.nrrd'

    annotation_file = os.path.join(annotation_path)
    if with_nrrd:
        annotation, _ = nrrd.read(annotation_file)
    else:
        # I'm not sure if anyones workflows use this so leaving it as an option, but 
        # making with_nrrd a required kwarg
        annotation = sitk.ReadImage( annotation_file )
    return annotation

def load_structure_graph():
    """
        Open up CCF structure graph data frame from disk

        typical protocol would be:
        cache = ReferenceSpaceCache(
        manifest=os.path.join("allen_ccf", "manifest.json"),  # downloaded files are stored relative to here
        resolution=10,
        reference_space_key="annotation/ccf_2017"  # use the latest version of the CCF
        )
        rsp = cache.get_reference_space()
        sg = rsp.remove_unassigned()
        sg_df = pd.DataFrame.from_records(sg)

    """
    sg_path =  files('morph_utils') / 'data/ccf_structure_graph.csv'
    df = pd.read_csv(sg_path)
    df['structure_id_path'] = df['structure_id_path'].apply(ast.literal_eval)
    df['structure_set_ids'] = df['structure_set_ids'].apply(ast.literal_eval)
    df['rgb_triplet'] = df['rgb_triplet'].apply(ast.literal_eval)
    df = df.set_index('acronym')
    return df


def process_pin_jblob( slide_specimen_id, jblob, annotation, structures, prints=False) :
    """
    Get CCF coordinates and structure for pins made with Cell Locator tool (starting mid 2022).

    :param slide_specimen_id: id of slide containing pins
    :param jblob: dictionary of pins for this slide made with the Cell Locator tool
    :param annotation: CCF annotation volume
    :param structures: DataFrame of all structures in CCF
    :return: list of dicts containing CCF location and structure of each pin in this slide
    """
    
    locs = []
    for m in jblob['markups'] :

        info = {}
        info['slide_specimen_id'] = slide_specimen_id
        info['specimen_name'] = m['name'].strip()
        try: info['specimen_id'] = int(get_id_by_name(info['specimen_name']))
        except: info['specimen_id'] = -1

        if m['markup']['type'] != 'Fiducial' :
            continue
            
        if 'controlPoints' not in m['markup'] :
            if prints: print(info)
            if prints: print("WARNING: no control point found, skipping")
            continue
            
        if m['markup']['controlPoints'] == None :
            if prints: print(info)
            if prints: print("WARNING: control point list empty, skipping")
            continue
            
        if len(m['markup']['controlPoints']) > 1 :
            if prints: print(info)
            if prints: print("WARNING: more than one control point, using the first")

        #
        # Cell Locator is LPS(RAI) while CCF is PIR(ASL)
        #
        pos = m['markup']['controlPoints'][0]['position']
        info['x'] =  1.0 * pos[1]
        info['y'] = -1.0 * pos[2]
        info['z'] = -1.0 * pos[0]
        
        if (info['x'] < 0 or info['x'] > 13190) or \
            (info['y'] < 0 or info['y'] > 7990) or \
            (info['z'] < 0 or info['z'] > 11390) :
            if prints: print(info)
            if prints: print("WARNING: ccf coordinates out of bounds")
            continue
        
        # Read structure ID from CCF
        point = (info['x'], info['y'], info['z'])
        
        # -- this simply divides cooordinates by resolution/spacing to get the pixel index
        pixel = annotation.TransformPhysicalPointToIndex(point)
        sid = annotation.GetPixel(pixel)
        info['structure_id'] = sid
        
        if sid not in structures.index :
            if prints: print(info)
            if prints: print("WARNING: not a valid structure - skipping")
            continue
        
        info['structure_acronym'] = structures.loc[sid]['acronym']

        locs.append(info)

    return locs

def get_soma_structure_and_ccf_coords():
    """
    Get CCF location and structure of all pins (somas and fiducials) 
    made with Cell Locator tool (starting mid 2022).

    :return: DataFrame containing CCF x,y,z coords and structure for all pins 
    """

    # (1) Get structure information from LIMS - this is only needed for validataion
    structures = get_structures()
    structures = pd.DataFrame.from_dict(structures)
    structures.set_index('id', inplace=True)

    # (2) Open up CCF annotation volume
    annotation = open_ccf_annotation(with_nrrd=False)

    # (3) Get json blobs (pin info) for all slides that have pins with Cell Locator tool
    pins = query_pinning_info_cell_locator()
    pins = pd.DataFrame.from_dict(pins)

    # (4) For each cell, convert Cell Locator to CCF coordinates and find structure using CCF annotation
    cell_info = []
    for index, row in pins.iterrows() :    
        jblob = row['data']
        processed = process_pin_jblob( row['specimen_id'], jblob, annotation, structures )
        cell_info.extend(processed)
    # (5) Return output as DataFrame
    df = pd.DataFrame(cell_info)
    return df

def move_soma_to_left_hemisphere(morph, resolution, volume_shape, z_midline):
    """
    Move a ccf registered morphology to the left hemisphere.

    Args:
        morph (Morphology): input morphology object (neuron_morphology.Morphology)
        resolution (int): number of um per voxel
        volume_shape (tuple): shape of ccf atlas in voxels
        z_midline (int): micron location of z-midline

    Returns:
        Morphology: translated morphology object
    """
    z_size = volume_shape[2]*resolution
    original_morph = morph.clone()
    soma = morph.get_soma()
    soma_z = soma['z'] 
    if soma_z > z_midline:
        new_soma_z = int(z_size - soma_z)

        # center on it's soma
        to_origin = aff.from_list([1, 0, 0, 0, 1, 0, 0, 0, 1, -soma['x'], -soma['y'], -soma['z']])
        to_origin.transform_morphology(morph)

        # mirror in z
        z_mirror = aff.from_list([1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0])
        z_mirror.transform_morphology(morph)

        # move back to original x and y and out to new z
        to_new_location = aff.from_list(
            [1, 0, 0, 0, 1, 0, 0, 0, 1, int(original_morph.get_soma()['x']), int(original_morph.get_soma()['y']), new_soma_z])
        to_new_location.transform_morphology(morph)

    return morph

def coordinates_to_voxels(coords, resolution=(10, 10, 10)):
    """ Find the voxel coordinates of spatial coordinates

    Parameters
    ----------
    coords : array
        (n, m) coordinate array. m must match the length of `resolution`
    resolution : tuple, default (10, 10, 10)
        Size of voxels in each dimension

    Returns
    -------
    voxels : array
        Integer voxel coordinates corresponding to `coords`
    """

    if len(resolution) != coords.shape[1]:
        raise ValueError(
            f"second dimension of `coords` must match length of `resolution`; "
            f"{len(resolution)} != {coords.shape[1]}")

    if not np.issubdtype(coords.dtype, np.number):
        raise ValueError(f"coords must have a numeric dtype (dtype is '{coords.dtype}')")

    voxels = np.floor(coords / resolution).astype(int)
    return voxels

def get_ccf_structure(voxel, name_map=None, annotation=None, coordinate_to_voxel_flag=True):
    """ 
    Will return the structure name for a given voxel. If it is out of cortex, returns Out Of Cortex


    Args:
        voxel (list): voxel location
        name_map (dict): dictionary that maps ccf structure id to structure name
        annotation (array): 3 dimensional ccf annotation array.
        coordinate_to_voxel_flag (bool, optional): _description_. Defaults to True.
    """
    if annotation is None:
        annotation = open_ccf_annotation(with_nrrd=True)
    
    if name_map is None:
        name_map = NAME_MAP
            
    if coordinate_to_voxel_flag:
        voxel = coordinates_to_voxels(voxel.reshape(1, 3))[0]

    voxel = voxel.astype(int)
    volume_shape = (1320, 800, 1140)
    for dim in [0,1,2]:
        if voxel[dim] == volume_shape[dim]:
            voxel[dim] = voxel[dim]-1

        if voxel[dim] >= volume_shape[dim]:
            # print("Dimension {} was provided values {} that exceeds volume size {}".format(dim,voxel[dim], volume_shape))
            return "Out Of Cortex"

    structure_id = annotation[voxel[0], voxel[1], voxel[2]]
    if structure_id == 0:
        return "Out Of Cortex"
    
    return name_map[structure_id]

def projection_matrix_for_swc(input_swc_file, branch_count, annotation=None, 
                              annotation_path = None, volume_shape=(1320, 800, 1140),
                              resolution=10):
    """
    Given a swc file, quantify the projection matrix. That is the amount of axon in each structure. This function assumes
    there is equivalent internode spacing (i.e. the input swc file should be resampled prior to running this code). 


    Args:
        input_swc_file (str): path to swc file
        branch_count (bool): if True, will count number of branches instead of the number of axon nodes
        annotation (array, optional): 3 dimensional ccf annotation array. Defaults to None.
        annotation_path (str, optional): path to nrrd file to use (optional). Defaults to None.
        volume_shape (tuple, optional): the size in voxels of the ccf atlas (annotation volume). Defaults to (1320, 800, 1140).
        resolution (int, optional): resolution (um/pixel) of the annotation volume
        
    Returns:
        filename (str)
        
        specimen_projection_summary (dict): keys are strings of structures and values are the quantitiave projection
        values. Either axon length, or number numbe of nodes depending on branch_count.
        
        specimen_projection_summary_branch_and_tip (dict): keys are structures and values are the quantitiave projection
        values. Either axon length, or number numbe of nodes. This dict only returns
        structures where there is a branch or a tip node in that structure.

    """

    if annotation is None:
        if isinstance(annotation_path, str):
            if not os.path.exists(annotation_path):
                resolution = 10
                volume_shape=(1320, 800, 1140)
                print(f"WARNING: Annotation path provided does not exist, defaulting to 10um resolution, (1320,800, 1140) ccf.\n{annotation_path}")
                annotation_path = None
        annotation = open_ccf_annotation(with_nrrd=True, annotation_path=annotation_path)
        

    sg_df = load_structure_graph()
    name_map = NAME_MAP
    full_name_to_abbrev_dict = dict(zip(sg_df.name, sg_df.index))
    full_name_to_abbrev_dict['Out Of Cortex'] = 'Out Of Cortex'
    fiber_tracts_id = sg_df[sg_df['name'] == 'fiber tracts']['id'].iloc[0]
    fiber_tract_acronyms = sg_df[sg_df['structure_id_path'].apply(lambda x: fiber_tracts_id in x)].index

    ventricular_system_id = sg_df[sg_df['name'] == 'ventricular systems']['id'].iloc[0]
    vs_acronyms = sg_df[sg_df['structure_id_path'].apply(lambda x: ventricular_system_id in x)].index

    z_size = resolution * volume_shape[2]
    z_midline = z_size / 2

    morph = morphology_from_swc(input_swc_file)
    morph = move_soma_to_left_hemisphere(morph, resolution, volume_shape, z_midline)    
    spacing = get_node_spacing(morph)[0]

    nodes_to_annotate = [n for n in morph.nodes() if (n['type'] == 2)]
    # print("Nodes to annotate before branch filter:")
    # print(len(nodes_to_annotate))
    if branch_count:
        nodes_to_annotate = [n for n in nodes_to_annotate if len(morph.get_children(n)) > 1]
        spacing = 1

    # print("Nodes to annotate:")
    # print(len(nodes_to_annotate))
    coords_to_annotate = np.array([[n['x'], n['y'], n['z']] for n in nodes_to_annotate])

    nodes_to_annotate_dict = {tuple([n['x'], n['y'], n['z']]): n['id'] for n in nodes_to_annotate}

    ipsi_coords = coords_to_annotate[coords_to_annotate[:, 2] < z_midline]
    contra_coords = coords_to_annotate[coords_to_annotate[:, 2] > z_midline]
    
    prefixes = {"ipsi": ipsi_coords,
                "contra": contra_coords}

    specimen_projection_summary = {}
    specimen_projection_summary_branch_and_tip = {}
    for prefix, coords_arr in prefixes.items():

        these_nodes = [morph.node_by_id(nodes_to_annotate_dict[tuple(c)]) for c in coords_arr]
        tip_and_branch_mask = [False] * len(these_nodes)
        for ct, no in enumerate(these_nodes):
            if len(morph.get_children(no)) != 1:
                tip_and_branch_mask[ct] = True

        # For each coordinate, get the ccf structure (full name with layer), abbreviate it
        structures = [full_name_to_abbrev_dict[get_ccf_structure(c, name_map, annotation, True)] for c in coords_arr]
        # add prefix and de-layer projection targets
        structures = [prefix + "_" + s for s in structures]
        projection_target_counts = pd.Series(structures).value_counts().to_dict()

        branch_and_tip_structures = list(set(np.array(structures)[tip_and_branch_mask]))
        # so that the nomenclature agrees with all projections
        branch_and_tip_structures = [s for s in branch_and_tip_structures]

        # Sort out fiber tracts
        curr_keys = list(projection_target_counts.keys())
        for projection_target in curr_keys:

            projection_value = projection_target_counts[projection_target]
            acronym = projection_target.replace(f"{prefix}_", "")
            if acronym in fiber_tract_acronyms:
                fiber_tract_key = f"{prefix}_fiber tracts"

                branch_and_tip_structures = list(
                    map(lambda x: x.replace(projection_target, fiber_tract_key), branch_and_tip_structures))

                if fiber_tract_key not in list(projection_target_counts.keys()):
                    projection_target_counts[fiber_tract_key] = 0
                projection_target_counts[fiber_tract_key] += projection_value

                del projection_target_counts[projection_target]

        ventral_targs = ["{}_{}".format(prefix, v) for v in vs_acronyms]
        targets_to_remove = [f"{prefix}_Out Of Cortex", f"{prefix}_root"] + ventral_targs
        for targ in targets_to_remove:
            if targ in projection_target_counts.keys():
                del projection_target_counts[targ]

            if targ in branch_and_tip_structures:
                branch_and_tip_structures.remove(targ)

        # Add them to bilateraldict
        for k, v in projection_target_counts.items():
            specimen_projection_summary[k] = v * spacing

            if k in branch_and_tip_structures:
                specimen_projection_summary_branch_and_tip[k] = v * spacing

    return input_swc_file, specimen_projection_summary, specimen_projection_summary_branch_and_tip