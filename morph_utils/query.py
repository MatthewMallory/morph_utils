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
