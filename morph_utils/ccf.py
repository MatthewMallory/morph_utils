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
import warnings
from copy import copy
import matplotlib.pyplot as plt
from morph_utils.query import get_id_by_name, get_structures, query_pinning_info_cell_locator
from morph_utils.measurements import get_node_spacing
from morph_utils.modifications import resample_morphology


NAME_MAP_FILE = files('morph_utils') / 'data/ccf_structure_name_map.json'
with open(NAME_MAP_FILE, "r") as fn: 
    NAME_MAP = json.load(fn)
NAME_MAP = {int(k):v for k,v in NAME_MAP.items()}

ACR_MAP_FILE = files('morph_utils') / 'data/ccf_structure_acronym_by_id.json'
with open(ACR_MAP_FILE, "r") as fn: 
    ACRONYM_MAP = json.load(fn)
ACRONYM_MAP = {k:int(v) for k,v in ACRONYM_MAP.items()}


_cached_ccf_annotation = None

def get_cached_ccf_annotation(annotation_file):
    """Retrieve the cached CCF annotation or load it if not already cached"""
    global _cached_ccf_annotation
    default_annotation_path = files('morph_utils') / 'data/annotation_10.nrrd'
    if (_cached_ccf_annotation is not None) and (annotation_file==default_annotation_path):
        return _cached_ccf_annotation
    else:
        _cached_ccf_annotation, _ = nrrd.read(annotation_file)
        
    return _cached_ccf_annotation


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
        annotation = get_cached_ccf_annotation(annotation_path)
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

def projection_matrix_for_swc(input_swc_file, mask_method = "tip_and_branch", 
                              tip_count = False, annotation=None, 
                              annotation_path = None, volume_shape=(1320, 800, 1140),
                              resolution=10, node_type_list=[2],
                              resample_spacing=None):
    """
    Given a swc file, quantify the projection matrix. That is the amount of axon in each structure. This function assumes
    there is equivalent internode spacing (i.e. the input swc file should be resampled prior to running this code). 


    Args:
        input_swc_file (str): path to swc file
        mask_method (str): method used to mask structures. If 'None', will return a projection matrix of all nodes. If 
        'tip_and_branch' will return a projection matrix masking only structures with tip and branch nodes. If 'tip'
        will only look at structures with tip nodes. And last, if 'branch' will only look at structures with 
        branch nodes.
        tip_count (bool): if True, will count number of tips instead of the length of axon
        annotation (array, optional): 3 dimensional ccf annotation array. Defaults to None.
        annotation_path (str, optional): path to nrrd file to use (optional). Defaults to None.
        volume_shape (tuple, optional): the size in voxels of the ccf atlas (annotation volume). Defaults to (1320, 800, 1140).
        resolution (int, optional): resolution (um/pixel) of the annotation volume
        node_type_list (list of ints): node type to extract projection data for, typically axon (2)
        resample_spacing (float or None): if not None, will resample the input morphology to the designated 
        internode spacing
        
    Returns:
        filename (str)
        
        specimen_projection_summary (dict): keys are strings of structures and values are the quantitiave projection
        values. Either axon length, or number numbe of nodes depending on tip_count argument.

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
    if resample_spacing is not None:
        morph = resample_morphology(morph, resample_spacing)
               
    spacing = get_node_spacing(morph)[0]

    morph_df = pd.DataFrame(morph.nodes())

    # filter by axon/dend types
    morph_df = morph_df[morph_df['type'].isin(node_type_list)]
    
    # annotate each node
    if morph_df.empty:
        print("Its empty")
        return input_swc_file, {} 
    
    morph_df['ccf_structure'] = morph_df.apply(lambda rw: full_name_to_abbrev_dict[get_ccf_structure( np.array([rw.x, rw.y, rw.z]) , name_map, annotation, True)], axis=1)

    # roll up fiber tracts
    morph_df.loc[morph_df['ccf_structure'].isin(fiber_tract_acronyms),'ccf_structure']='fiber tracts'
    
    # identify branch/tip/reducible
    def node_ider(morph,i):
        nc = len(morph.child_ids([i])[0]) 
        if nc==0:
            return 'tip'
        elif nc>1:
            return 'branch'
        else:
            return 'reducible'
    morph_df["node_type"] = morph_df.id.apply(lambda i: node_ider(morph,i))

    # determine ipsi/contra projections
    morph_df["ccf_structure_sided"] = morph_df.apply(lambda row: "ipsi_{}".format(row.ccf_structure) if row.z<z_midline else "contra_{}".format(row.ccf_structure), axis=1)


    # mask the morphology dataframe accordinagly
    if mask_method!="None":
            
        keep_structs = []
        for struct, struct_df in morph_df.groupby("ccf_structure_sided"):
            node_types_in_struct = struct_df.node_type.unique().tolist()
        
            if (mask_method == 'tip') & ("tip" in node_types_in_struct):
                keep_structs.append(struct)
                
            elif (mask_method == 'branch') & ("branch" in node_types_in_struct):
                keep_structs.append(struct)
                
            elif (mask_method == 'tip_and_branch') & (all([i in node_types_in_struct for i in ['tip','branch']])):
                keep_structs.append(struct)
                
            elif (mask_method == 'tip_or_branch') & (any([i in node_types_in_struct for i in ['tip','branch']])):
                keep_structs.append(struct)
                            
        morph_df_masked = morph_df[morph_df['ccf_structure_sided'].isin(keep_structs)]
        
    else:
        morph_df_masked = morph_df

    # remove ventral targets and out of brain 
    ventral_targs = ["ipsi_{}".format(v) for v in vs_acronyms] + ["contra_{}".format(v) for v in vs_acronyms]
    targets_to_remove = ["ipsi_Out Of Cortex", "ipsi_root","contra_Out Of Cortex", "contra_root"] + ventral_targs
    morph_df_masked = morph_df_masked[~morph_df_masked['ccf_structure_sided'].isin(targets_to_remove)]

    # accomodate tip counting instead of axon length 
    if tip_count:
        morph_df_masked = morph_df_masked[morph_df_masked['node_type']=='tip']
        spacing = 1
    
    # qunatify projections per structure 
    n_nodes_per_structure = morph_df_masked.ccf_structure_sided.value_counts()
    axon_length_per_structure = n_nodes_per_structure*spacing
    specimen_projection_summary = axon_length_per_structure.to_dict()
    
    return input_swc_file, specimen_projection_summary   

 
def correct_superficial_nodes_out_of_brain(morphology,
                                           annotation,
                                           closest_surface_voxel_file,
                                           surface_paths_file,
                                           tree,
                                           volume_shape=(1320, 800, 1140),
                                           resolution=10,
                                           isocortex_struct_id=315,
                                           generate_plot=True,
                                           fig_ofile=None,
                                           ):
    """
    This function attempts to correct nodes that appear out of brain due to registrastion
    issues. This will find the streamline that passes closest to the cells soma and 
    slide the cell depper along that streamline until stopping conditions have been satisifed.
    Where stopping conditions are either all the nodes are in the brain, or the cell has been
    pushed to the bottom of the streamline. 
    
    NOTE: this function should only be used on local morphologies. It is not recommended to apply this
    function to the entire cell. This function is on attempting to fix local issues for more accurate 
    local feature calcualtion. Local morphologies can be generated from:
    skeleton_keys.full_morph.local_crop_cortical_morphology 
    or 
    morph_utils.executable_scripts.local_crop_ccf_swc_directory)
    
    Args:
        morphology (neuron_morphology.Morphology): A LOCAL morphology (derived from full_morph.local_crop_cortical_morphology or morph_utils.executable_scripts.local_crop_ccf_swc_directory)
        annotation (3d np.array): ccf annotation atlas
        closest_surface_voxel_file (str): path to closest_surface_voxel_file
        surface_paths_file (str):  path to surface_paths_file
        tree (_type_): allensdk reference space tree
        volume_shape (tuple, optional): shape of annotation. Defaults to (1320, 800, 1140).
        resolution (int, optional): resolution of atlas. Defaults to 10.
        isocortex_struct_id (int, optional): structure id for isocortex. Defaults to 315.
        generate_plot (bool, optional): whether to generate qc plots or not. Defaults to True.
        fig_ofile (str, optional): path to save qc plot at. Defaults to None.
        
    Returns:
        tuple (morphology (neuron_morphology.Morphology), move_bool) return morphology and bool
        indicating if the morphology was moved
    """
    
    
    from sklearn.neighbors import KDTree
    try:   
        from skeleton_keys import full_morph
    except ImportError:
        msg = """
        Required module (skeleton_keys.full_morph) is not installed. It's possible you have skeleton_keys installed
        but not the correct branch/version. As of 12/22/23 the full_morph branch has not been merged into the main 
        branch of skeleton_keys so check the full_morph-MM-edits branch for the full_morph features used in this code.
        This module is only needed for the function morph_utils.ccf.correct_superficial_nodes_out_of_brain. If you are
        not using this function, no need to install skeleton-keys.
        """
        warnings.warn(msg)
    try:
        from ccf_streamlines.angle import find_closest_streamline
    except ImportError:
        msg = """
        ccf_streamlines is required for this function. Please reference the link below for installation.
        
        https://github.com/AllenInstitute/ccf_streamlines
        """
        warnings.warn(msg)

    
    
    
    morph = morphology.clone()
    
    morph_coords = np.array([ [n['x'], n['y'], n['z'] ] for n in morph.nodes()])
    morph_voxels = coordinates_to_voxels(morph_coords)
    out_of_brain_voxels = [v for v in morph_voxels if annotation[v[0],v[1],v[2]] == 0]  
    
    if not out_of_brain_voxels:
        return morph,False
    else:
        
        # find the streamline closest to the soma
        # if the cells soma is in WM like some deep L6bs, 
        # we cannot push the cell any deeper. This approach uses streamlines
        # and streamlines do not extend into WM so we do not have an orientation on
        # how to push those cells. 
        cells_soma = morph.get_soma()
        soma_arr = np.array([cells_soma['x'],cells_soma['y'], cells_soma['z']]).reshape(1,3)
        soma_out_of_cortex_bool, nearest_cortex_coord = full_morph.check_coord_out_of_cortex(soma_arr,
                                                                                    isocortex_struct_id,
                                                                                    atlas_volume=annotation,
                                                                                    closest_surface_voxel_file=closest_surface_voxel_file,
                                                                                    surface_paths_file=surface_paths_file,
                                                                                    tree=tree)

        if soma_out_of_cortex_bool:
            msg = """WARNING: Can not correct out of brain nodes. Unable to identify streamline
            nearest to the cells soma because the soma is located out of cortex (likely in white matter)
            """
            warnings.warn(msg)
            return morph
        
        # original_soma_arr = copy(soma_arr)
        closest_streamline = find_closest_streamline(soma_arr,
                            closest_surface_voxel_file,
                            surface_paths_file,
                            resolution=(10,10,10),
                            volume_shape=volume_shape
                           )
        
        streamline_kd_tree = KDTree(closest_streamline)

        # find streamline node closest to the soma
        dist, streamline_indices = streamline_kd_tree.query(soma_arr)
        streamline_index = streamline_indices[0][0]
        nearest_streamline_node = closest_streamline[streamline_index]

        # this transform is what we will apply every step we move down the streamline
        # so we will move one node down the streamline -> apply this transofrm -> check out of brain nodes -> repeat
        deltas_from_streamline = soma_arr[0] - nearest_streamline_node
        dx, dy, dz = deltas_from_streamline[0], deltas_from_streamline[1], deltas_from_streamline[2]
        aff_from_streamline = [1, 0, 0, 0, 1, 0, 0, 0, 1, dx, dy, dz]
        offset_transformation = aff.from_list(aff_from_streamline)


        # positive 1 to move down/deeper along the streamline
        index_mover = 1

        stopping_condition = False
        while stopping_condition == False:
            
            # move one streamline 
            streamline_index += index_mover
            if streamline_index >= len(closest_streamline)-1:
                stopping_condition=True
                warn_msg = """
                WARNING, cell has been moved to the end of the streamline, but there are still {}
                nodes out of brain.""".format(len(out_of_brain_voxels))
                warnings.warn(warn_msg)
                break
                
            current_streamline_node_to_check = closest_streamline[streamline_index]

            # move cell to next streamline
            deltas_to_current_streamline = current_streamline_node_to_check - soma_arr[0]
            dx_curr, dy_curr, dz_curr = deltas_to_current_streamline[0], deltas_to_current_streamline[1], deltas_to_current_streamline[2]
            
            aff_to_current_streamline = [1, 0, 0, 0, 1, 0, 0, 0, 1, dx_curr, dy_curr, dz_curr]

            # apply offset
            aff.from_list(aff_to_current_streamline).transform_morphology(morph)
            offset_transformation.transform_morphology(morph)

            # update soma 
            this_soma = morph.get_soma()
            soma_arr = np.array([this_soma['x'], this_soma['y'], this_soma['z']]).reshape(1,3)

            # measure coordinates that are still out of brain
            this_morph_coords = np.array([ [n['x'], n['y'], n['z'] ] for n in morph.nodes()])
            this_morph_voxels = coordinates_to_voxels(this_morph_coords)
            out_of_brain_voxels = [v for v in this_morph_voxels if annotation[v[0],v[1],v[2]] == 0]

            if out_of_brain_voxels == []:
                stopping_condition = True
                
    if generate_plot:
                
        streamline_vox = coordinates_to_voxels(closest_streamline)
        fig,axes=plt.subplots(3,1)
        crops=[False,True]
        for axe,crop in zip(axes[:-1],crops):
            
            soma_x = soma_arr[0][0]*(1/resolution)
            soma_y = soma_arr[0][1]*(1/resolution)
            soma_z = soma_arr[0][2]*(1/resolution)
            
            atlas_slice = annotation[int(soma_x),:,:].astype(bool)
            axe.imshow(atlas_slice)
            axe.scatter(this_morph_voxels[:,2],this_morph_voxels[:,1],s=0.1)
            axe.scatter(streamline_vox[:,2],streamline_vox[:,1],s=0.1)
            axe.scatter(soma_z,soma_y,s=10,c='r',marker='X')
            if crop:
                    
                buff=100
                axe.set_xlim(soma_z-buff,soma_z+buff)
                axe.set_ylim(soma_y-buff,soma_y+buff)

        axe = axes[2]
        axe.scatter(morph_voxels[:,2],this_morph_voxels[:,1],s=0.5,alpha=0.75,label='original morph')
        axe.scatter(this_morph_voxels[:,2],this_morph_voxels[:,1],s=0.5,alpha=0.75,label='moved morph')
        axe.plot(streamline_vox[:,2],streamline_vox[:,1],lw=3,c='g')
        axe.legend()
        axe.set_aspect('equal')
        fig.set_size_inches(5,12)   
        if fig_ofile is not None:
            fig.savefig(fig_ofile,dpi=300,bbox_inches='tight')
        plt.clf()
       
    return morph, True
