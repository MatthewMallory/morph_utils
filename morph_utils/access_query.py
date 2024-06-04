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
    result = query_access(db_file, query)
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
    result = query_access(db_file, query)
    return result


def get_nhp_bg_passing(db_file):
    """
    Get NHP cells in basal ganglia that can be reconstructed.

    :return: passing nhp basal ganglia cells 
    """
    query = """
       SELECT 
              [Cell Specimen Id],
              [Cell Overall State], 
              Project, 
              [Pinned Structure and Layer], 
              NHP_BG_Mapping, 
              NHP_VIS_Mapping, 
              NHP_MTG_Mapping
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              (Project LIKE 'qIVSCC-MET%' OR Project = 'MET-NM')
              AND ([Cell Overall State] NOT LIKE 'Deferred' 
              AND [Cell Overall State] NOT LIKE 'To be%' 
              AND [Cell Overall State] NOT LIKE 'QC%' 
              AND [Cell Overall State] NOT LIKE 'Failed%' 
              AND [Cell Overall State] NOT LIKE 'Rescan%')
              AND ([Pinned Structure and Layer] LIKE 'STr%' 
                     OR [Pinned Structure and Layer] LIKE 'dSTR%' 
                     OR [Pinned Structure and Layer] LIKE 'Ca%' 
                     OR [Pinned Structure and Layer] LIKE 'Pu%' 
                     OR [Pinned Structure and Layer] LIKE 'vSTr%' 
                     OR [Pinned Structure and Layer] LIKE 'NAC%' 
                     OR [Pinned Structure and Layer] LIKE 'ISC%' 
                     OR [Pinned Structure and Layer] LIKE 'Tu%' 
                     OR [Pinned Structure and Layer] LIKE 'GP%' 
                     OR [Pinned Structure and Layer] LIKE 'IC%');
    """
    result = query_access(db_file, query)
    return result


def get_nhp_cortical_passing(db_file):
    """
    Get NHP cortical cells that can be reconstructed.

    :return: passing nhp cortical cells 
    """
    query = """
       SELECT 
              [Cell Specimen Id], 
              [Cell Overall State],
              Project, 
              [Pinned Structure and Layer], 
              NHP_BG_Mapping, 
              NHP_VIS_Mapping, 
              NHP_MTG_Mapping
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              (Project LIKE 'qIVSCC-MET%' OR Project = 'MET-NM')
              AND ([Cell Overall State] NOT LIKE 'Deferred' 
              AND [Cell Overall State] NOT LIKE 'To be%' 
              AND [Cell Overall State] NOT LIKE 'QC%' 
              AND [Cell Overall State] NOT LIKE 'Failed%' 
              AND [Cell Overall State] NOT LIKE 'Rescan%')
              AND ([Pinned Structure and Layer] LIKE 'NC%' 
                     OR [Pinned Structure and Layer] LIKE 'dlPF%' 
                     OR [Pinned Structure and Layer] LIKE 'mOF%' 
                     OR [Pinned Structure and Layer] LIKE 'cOF%' 
                     OR [Pinned Structure and Layer] LIKE 'MFC%' 
                     OR [Pinned Structure and Layer] LIKE 'rCG%' 
                     OR [Pinned Structure and Layer] LIKE 'Oc%' 
                     OR [Pinned Structure and Layer] LIKE 'V1%' 
                     OR [Pinned Structure and Layer] LIKE 'V2%' 
                     OR [Pinned Structure and Layer] LIKE 'HIP%'
                     OR [Pinned Structure and Layer] LIKE 'M1C%'
                     OR [Pinned Structure and Layer] LIKE 'STG%');
    """
    result = query_access(db_file, query)
    return result


def get_nhp_other_passing(db_file):
     """
     Get NHP cells from regions outside of basal ganglia and cortex that can be reconstructed.
     
     :return: passing nhp cells from regions outside of basal ganglia and cortex  
     """
     passing = get_nhp_passing(db_file)
     bg_passing = get_nhp_bg_passing(db_file)
     cortex_passing = get_nhp_cortical_passing(db_file)
     other_passing = passing[~passing['Cell Specimen Id'].isin(bg_passing['Cell Specimen Id']) & ~passing['Cell Specimen Id'].isin(cortex_passing['Cell Specimen Id'])]
     
     return other_passing


def get_mouse_passing(db_file):
    """
    Get mouse cells that can be reconstructed. 

    :return: passing mouse cells 
    """

    query = """
    SELECT 
       [Cell Specimen Id], 
       [Cell Overall State], 
       Project, 
       [Pinned Structure and Layer], 
       mouse_wholebrain_mapping,
       mouse_wholebrain_supertype
    FROM IVSCCTrackingDatabaseProduction
    WHERE 
       (
              Project LIKE 'mIVSCC-MET' OR 
              Project LIKE 'T301%'
       ) AND 
       (
              [Cell Overall State] NOT LIKE 'Deferred' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%'
       );
    """
    
    result = query_access(db_file, query)
    return result


def get_mouse_bg_passing(db_file):
    """
    Get mouse cells in basal ganglia that can be reconstructed. 

    :return: passing mouse cells in basal ganglia 
    """
    query = """ 
       SELECT 
              [Cell Specimen Id], 
              [Cell Overall State],
              Project, 
              [Pinned Structure and Layer], 
              mouse_wholebrain_mapping,
              mouse_wholebrain_supertype
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
              (
                     Project LIKE 'mIVSCC-MET' OR 
                     Project LIKE 'T301%'
              ) AND 
              (
                     [Cell Overall State] NOT LIKE 'Deferred' AND 
                     [Cell Overall State] NOT LIKE 'To be%' AND 
                     [Cell Overall State] NOT LIKE 'QC' AND 
                     [Cell Overall State] NOT LIKE 'Failed%' AND 
                     [Cell Overall State] NOT LIKE 'Rescan%'
              ) AND 
              (
                     [Pinned Structure and Layer] LIKE 'STR%' OR 
                     [Pinned Structure and Layer] = 'CP' OR 
                     [Pinned Structure and Layer] = 'ACB' OR 
                     [Pinned Structure and Layer] = 'FS' OR 
                     [Pinned Structure and Layer] LIKE 'OT%' OR 
                     [Pinned Structure and Layer] LIKE 'isl%' OR 
                     [Pinned Structure and Layer] = 'islm' OR 
                     [Pinned Structure and Layer] = 'LSS' OR 
                     [Pinned Structure and Layer] = 'PALd' OR 
                     [Pinned Structure and Layer] = 'GPe' OR 
                     [Pinned Structure and Layer] = 'GPi' OR 
                     [Pinned Structure and Layer] = 'SN' OR 
                     [Pinned Structure and Layer] = 'SNc' OR 
                     [Pinned Structure and Layer] = 'SNr' OR 
                     [Pinned Structure and Layer] = 'STN' OR 
                     [Pinned Structure and Layer] LIKE 'PST%' OR 
                     [Pinned Structure and Layer] LIKE 'PSTN%'
              );

    """
    result = query_access(db_file, query)
    return result


def get_mouse_hpf_passing(db_file):
    """
    Get mouse cells in hippocampal formation that can be reconstructed. 

    :return: passing mouse cells in hippocampal formation  
    """
    query = """
       SELECT 
              [Cell Specimen Id], 
              [Cell Overall State], 
              Project, 
              [Pinned Structure and Layer], 
              mouse_wholebrain_mapping,
              mouse_wholebrain_supertype
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
       (
              Project LIKE 'mIVSCC-MET' 
       ) AND 
       (
              [Cell Overall State] NOT LIKE 'Deferred' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%'
       ) AND 
       (
              [Pinned Structure and Layer] LIKE 'HPF%' OR 
              [Pinned Structure and Layer] LIKE 'HIP%' OR 
              [Pinned Structure and Layer] LIKE 'CA%' OR 
              [Pinned Structure and Layer] LIKE 'DG%' OR 
              [Pinned Structure and Layer] LIKE 'FC%' OR 
              [Pinned Structure and Layer] LIKE 'IG%' OR 
              [Pinned Structure and Layer] LIKE 'RHP%' OR 
              [Pinned Structure and Layer] LIKE 'ENT%' OR 
              [Pinned Structure and Layer] LIKE 'PAR%' OR 
              [Pinned Structure and Layer] LIKE 'POST%' OR 
              [Pinned Structure and Layer] LIKE 'PRE%' OR 
              [Pinned Structure and Layer] LIKE 'SUB%' OR 
              [Pinned Structure and Layer] LIKE 'ProS%' OR 
              [Pinned Structure and Layer] LIKE 'HATA%' OR 
              [Pinned Structure and Layer] LIKE 'APr%'
       );
    """
    result = query_access(db_file, query)
    return result


def get_mouse_isocortex_passing(db_file):
    """
    Get mouse cells in isocortex that can be reconstructed. 

    :return: passing mouse cells in isocortex
    """
    query = """
       SELECT               
              [Cell Specimen Id], 
              [Cell Overall State], 
              Project, 
              [Pinned Structure and Layer], 
              mouse_wholebrain_mapping,
              mouse_wholebrain_supertype
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
       (
              Project LIKE 'mIVSCC-MET' OR 
              Project LIKE 'T301%'
       ) AND 
       (
              [Cell Overall State] NOT LIKE 'Deferred' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%'
       ) AND 
       (
              [Pinned Structure and Layer] LIKE 'FRP%' OR 
              [Pinned Structure and Layer] LIKE 'MO%' OR 
              [Pinned Structure and Layer] LIKE 'SS%' OR 
              [Pinned Structure and Layer] LIKE 'GU%' OR 
              [Pinned Structure and Layer] LIKE 'VISC%' OR 
              [Pinned Structure and Layer] LIKE 'AUD%' OR 
              [Pinned Structure and Layer] LIKE 'VIS%' OR 
              [Pinned Structure and Layer] LIKE 'ACA%' OR 
              [Pinned Structure and Layer] LIKE 'PL%' OR 
              [Pinned Structure and Layer] LIKE 'ILA%' OR 
              [Pinned Structure and Layer] LIKE 'ORB%' OR 
              [Pinned Structure and Layer] LIKE 'AI%' OR 
              [Pinned Structure and Layer] LIKE 'RSP%' OR 
              [Pinned Structure and Layer] LIKE 'PTLp%' OR 
              [Pinned Structure and Layer] LIKE 'TEa%' OR 
              [Pinned Structure and Layer] LIKE 'PERI%' OR 
              [Pinned Structure and Layer] LIKE 'ECT%' OR 
              [Pinned Structure and Layer] LIKE 'CTXsp%'
       );
    """
    result = query_access(db_file, query)
    return result


def get_mouse_thalamus_passing(db_file):
    """
    Get mouse cells in thalamus that can be reconstructed. 

    :return: passing mouse cells in thalamus
    """
    query = """
       SELECT 
              [Cell Specimen Id], 
              [Cell Overall State], 
              Project,
              [Pinned Structure and Layer], 
              mouse_wholebrain_mapping,
              mouse_wholebrain_supertype
       FROM IVSCCTrackingDatabaseProduction
       WHERE 
       (
              Project LIKE 'mIVSCC-MET' OR 
              Project LIKE 'T301%'
       ) AND 
       (
              [Cell Overall State] NOT LIKE 'Deferred' AND 
              [Cell Overall State] NOT LIKE 'To be%' AND 
              [Cell Overall State] NOT LIKE 'QC' AND 
              [Cell Overall State] NOT LIKE 'Failed%' AND 
              [Cell Overall State] NOT LIKE 'Rescan%'
       ) AND 
       (
              [Pinned Structure and Layer] IN (
              'TH', 'DORsm', 'VENT', 'VAL', 'VM', 'VP', '', 'PoT', 
              'SPF', 'SPFm', 'SPFp', 'SPA', 'PP', 'GENd', 'LAT', 
              'LP', 'PO', 'POL', 'SGN', 'Eth', 'REth', 'ATN', 'AV', 
              'AD', 'IAM', 'IAD', 'LD', 'MED', 'IMD', 'SMT', 'PR', 
              'MTN', 'PVT', 'PT', 'RE', 'Xi', 'ILM', 'RH', 'CM', 
              'PCN', 'CL', 'PF', 'PIL', 'RT', 'GENv', 'IGL', 'IntG', 
              'SubG', 'EPI', 'MH', 'LH', 'PIN'
              )
              OR 
              [Pinned Structure and Layer] LIKE 'VPL%' OR 
              [Pinned Structure and Layer] LIKE 'VPM%' OR 
              [Pinned Structure and Layer] LIKE 'MG%' OR 
              [Pinned Structure and Layer] LIKE 'LGd%' OR 
              [Pinned Structure and Layer] LIKE 'AM%' OR 
              [Pinned Structure and Layer] LIKE 'MD%' OR 
              [Pinned Structure and Layer] LIKE 'LGv%'
       );
    """    
    result = query_access(db_file, query)
    return result


def get_mouse_other_passing(db_file):
     """
     Get mouse cells from regions outside of basal ganglia, hippocampal formation, isocortex, and thalamus that can be reconstructed.
     
     :return: passing mouse cells from regions outside of basal ganglia, hippocampal formation, isocortex, and thalamus  
     """
     passing = get_mouse_passing(db_file)
     bg_passing = get_mouse_bg_passing(db_file)
     hpf_passing = get_mouse_hpf_passing(db_file)
     isocortex_passing = get_mouse_isocortex_passing(db_file)
     thalamus_passing = get_mouse_thalamus_passing(db_file)
     other_passing = passing[~passing['Cell Specimen Id'].isin(bg_passing['Cell Specimen Id']) & 
                             ~passing['Cell Specimen Id'].isin(hpf_passing['Cell Specimen Id']) &
                             ~passing['Cell Specimen Id'].isin(isocortex_passing['Cell Specimen Id']) &
                             ~passing['Cell Specimen Id'].isin(thalamus_passing['Cell Specimen Id'])]
     
     return other_passing
