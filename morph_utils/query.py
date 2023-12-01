import os
import allensdk.internal.core.lims_utilities as lu
from functools import partial

def default_query_engine():
    """Get Postgres query engine with environmental variable parameters"""

    return partial(
        lu.query,
        host=os.getenv("LIMS_HOST"),
        port=5432,
        database=os.getenv("LIMS_DBNAME"),
        user=os.getenv("LIMS_USER"),
        password=os.getenv("LIMS_PASSWORD")
    )

def get_name_by_id(specimen_id, query_engine=None):
    """
    Get cell name from id

    :param specimen_id: cell specimen id
    :return: cell specimen name 
    """
    if query_engine is None: query_engine = default_query_engine()

    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.id = '{}'
    """.format(specimen_id)
    
    x = query_engine(sql)[0]['sp']
    return x

def get_id_by_name(specimen_name, query_engine=None):
    """
    Get cell id from name

    :param specimen_name: cell specimen name
    :return: cell specimen id 
    """
    if query_engine is None: query_engine = default_query_engine()

    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.name = '{}'
    """.format(specimen_name)
    
    x = query_engine(sql)[0]['id']
    return x

#convert px -> um 
def query_for_z_resolution(specimen_id, query_engine=None):
    """
    Get resolution of z slice axis for a cell. 

    :param specimen_id: cell specimen id
    :return: z axis resolution
    """
    if query_engine is None: query_engine = default_query_engine()
     
    sql = """
    select ss.id, ss.name, shs.thickness from specimens ss
    join specimen_blocks sb on ss.id = sb.specimen_id
    join blocks bs on bs.id = sb.block_id
    join thicknesses shs on shs.id = bs.thickness_id 
    where ss.id = {}
    """.format(specimen_id)
    
    res = query_engine(sql)
    try:
        return res[0]['thickness']
    except:
        return None

def get_structures(query_engine=None):
    """
    Get structure information (from LIMS)
    
    :return: list of dicts with info for each brain structure
    """
    if query_engine is None: 
        query_engine = default_query_engine()

    sql = """SELECT * FROM structures where ontology_id = 1"""
    structures = query_engine(sql)

    return structures

def query_pinning_info_cell_locator(query_engine=None):
    """
    Get CCF pins made with Cell Locator tool (starting mid 2022)
    
    :return: list of dicts containing pins for each slide
    """
    if query_engine is None: 
        query_engine = default_query_engine()

    sql = """SELECT sm.* FROM specimen_metadata sm WHERE sm.current = 't' AND sm.kind = 'IVSCC tissue review'"""
    pins = query_engine(sql)

    return pins


def query_pinning_info(project_codes=["T301", "T301x", "mIVSCC-MET"], query_engine=None):
    """
    Get the pinned CCF coordinates for a set of projects. 
    (pre switch to Cell Locator tool mid 2022) 
    """
    if query_engine is None:
        query_engine = default_query_engine()

    project_codes_str = ", ".join([f"'{s}'" for s in project_codes])
    query = f"""
        select distinct
            sp.id as specimen_id,
            csl.x as x,
            csl.y as y,
            csl.z as z,
            slice.id as slice_id,
            slab.id as slab_id,
            brain.id as brain_id,
            a3d.*
        from specimens sp
        join specimens slice on slice.id = sp.parent_id
        join specimens slab on slab.id = slice.parent_id
        join specimens brain on brain.id = slab.parent_id
        join alignment3ds a3d on slice.alignment3d_id = a3d.id
        join projects prj on prj.id = sp.project_id
        left join cell_soma_locations csl on sp.id = csl.specimen_id
        where prj.code in ({project_codes_str})
    """
    results = query_engine(query)
    return results

def get_swc_from_lims(specimen_id, query_engine=None):

    if query_engine is None:
        query_engine = default_query_engine()

    query = "SELECT f.filename, f.storage_directory FROM \
     neuron_reconstructions n JOIN well_known_files f ON n.id = f.attachable_id \
     AND n.specimen_id = {} AND n.manual AND NOT n.superseded AND f.well_known_file_type_id = 303941301".format(specimen_id)
    
    result = query_engine(query)

    if result is None:
        raise("No SWC file found for specimen ID {}".format(specimen_id))
        
    swc_filename = result[0]['filename']
    swc_path = result[0]['storage_directory'] + swc_filename

    return swc_filename, swc_path