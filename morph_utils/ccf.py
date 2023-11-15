import os
import pandas as pd
from importlib.resources import files
import SimpleITK as sitk
from morph_utils.query import get_id_by_name, get_structures, query_pinning_info_cell_locator

def open_ccf_annotation(annotation_path=None, as_array=False):
    """
        Open up CCF annotation volume
    """
    if annotation_path is None: 
        annotation_path =  files('morph_utils') / 'data/annotation_10.nrrd'

    annotation_file = os.path.join(annotation_path)
    annotation = sitk.ReadImage( annotation_file )
    if as_array:
        return sitk.GetArrayFromImage(annotation)
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
    annotation = open_ccf_annotation()

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

