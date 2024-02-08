import pyodbc
import pandas as pd


def query_access(db_file, query):
    """
    Query Access database.

    :return: query result as pandas dataframe 
    """
    conn = pyodbc.connect(r'Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;' % db_file)
    result = pd.read_sql(query, conn)
    conn.close()
    return result


def get_nhp_passing(db_file):
    """
    Get NHP cells that can be reconstructed.

    :return: passing nhp cells 
    """
    query = """
       SELECT [Cell Specimen Id],
              [Cell Overall State], 
              Project, 
              [Pinned Structure and Layer], 
              NHP_BG_Mapping, 
              NHP_VIS_Mapping, 
              NHP_MTG_Mapping
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              (Project LIKE 'qIVSCC-MET%' OR Project = 'MET-NM') AND
              [Cell Overall State] NOT LIKE 'Deferred' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC%' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%';
    """
    nhp_passing = query_access(db_file, query)
    return nhp_passing


def get_human_passing(db_file):
    """
    Get human cells that can be reconstructed. 
    **These will always be cortical cells**

    :return: passing human cells 
    """
    query="""
       SELECT 
              [Cell Specimen Id],
              [Cell Overall State], 
              Project, 
              [Pinned Structure and Layer], 
              [Tree Mapping]
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              Project LIKE 'h%' AND
              [Cell Overall State] NOT LIKE 'Deferred%' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC%' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%' AND 
              [Cell Overall State] NOT LIKE 'Reconstruction IP - 0% (truncated)' AND 
              [Cell Overall State] NOT LIKE 'Uploaded (truncated)';
    """
    human_passing = query_access(db_file, query)
    return human_passing


def get_mouse_subcortical_passing(db_file):
    """
    Get mouse subcortical cells that can be reconstructed. 

    :return: passing mouse subcortical cells 
    """
    query = """
    SELECT  
       [Cell Specimen Id], 
       [Cell Overall State], 
       Project,
       [Pinned Structure and Layer], 
       mouse_wholebrain_mapping
    FROM IVSCCTrackingDatabaseProduction
    WHERE 
       Project LIKE 'mI%' AND
       [Cell Overall State] NOT LIKE 'Deferred' AND 
       [Cell Overall State] NOT LIKE 'To be%' AND 
       [Cell Overall State] NOT LIKE 'QC%' AND 
       [Cell Overall State] NOT LIKE 'Failed%' AND 
       [Cell Overall State] NOT LIKE 'Rescan%' AND 
       (
              ([Pinned Structure and Layer] NOT LIKE '%1%' AND 
              [Pinned Structure and Layer] NOT LIKE '%2/3%' AND 
              [Pinned Structure and Layer] NOT LIKE '%4%' AND 
              [Pinned Structure and Layer] NOT LIKE '%5%' AND 
              [Pinned Structure and Layer] NOT LIKE '%6%' AND 
              [Pinned Structure and Layer] NOT LIKE 'Vis%' AND 
              [Pinned Structure and Layer] NOT LIKE 'RSP%' AND 
              [Pinned Structure and Layer] NOT LIKE 'ORB%' AND 
              [Pinned Structure and Layer] NOT LIKE 'TEa%' AND 
              [Pinned Structure and Layer] NOT LIKE 'ENT%' AND 
              [Pinned Structure and Layer] NOT LIKE 'Ctx%' AND 
              [Pinned Structure and Layer] NOT LIKE 'MO')
              OR
              [Pinned Structure and Layer] Is Null
       );
    """
    mouse_subcortical_passing = query_access(db_file, query)
    return mouse_subcortical_passing


def get_mouse_dendrite_only_passing(db_file):
    """
    Get mouse dendrite only cells that can be reconstructed. 

    :return: passing mouse dendrite only cells 
    """
    query = """
       SELECT 
              Project,
              [Cell Overall State],
              [Cell Specimen Id],
              [Pinned Structure and Layer],
              [Tree Mapping]
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              Project LIKE 'm%' AND
              [Cell Overall State] = 'Reconstructable (DO)';
    """
    mouse_dendrite_only_passing = query_access(db_file, query)
    return mouse_dendrite_only_passing
