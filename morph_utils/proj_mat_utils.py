import os
import numpy as np
import pandas as pd

from morph_utils.ccf import STRUCTURE_DESCENDANTS_ACRONYM


def de_layer(st):
    CTX_STRUCTS = STRUCTURE_DESCENDANTS_ACRONYM['CTX']
    sub_st = st.replace("ipsi_","").replace("contra","")
    if sub_st in CTX_STRUCTS:
            
        for l in ["1","2/3","4","5","6a","6b"]:
            st = st.replace(l,"")
            
        if "ENT" in st:
            for l in ["2", "3","6"]:
                st = st.replace(l,"")
            
        return st
    else:
        return st


def roll_up_proj_mat(infile, outfile):
    
    df = pd.read_csv(infile, index_col=0)
    df.index = df.index.map(os.path.basename)
    
    non_proj_cols = [f for f in df.columns if not any([i in f for i in ["ipsi","contra"]])]
    new_df = df[non_proj_cols].copy()
    
    proj_cols = [f for f in df.columns if any([i in f for i in ["ipsi","contra"]])]
    de_layer_dict = {p:de_layer(p) for p in proj_cols}
    
    parent_names = list(de_layer_dict.values())
    unique_parent_names = np.unique(parent_names)
    unique_parent_names = sorted(unique_parent_names, key=lambda x:parent_names.index(x))
    
    roll_up_records = {}
    for low_res_struct in unique_parent_names:
        children = [k for k,v in de_layer_dict.items() if v==low_res_struct ]
        roll_up_records[low_res_struct] = children
    
    
    
    # for parent, child_list in roll_up_records.items():
    #     new_df[parent] = df[child_list].sum(axis=1)
    new_cols = {
        parent: df[child_list].sum(axis=1)
        for parent, child_list in roll_up_records.items()
    }
    new_cols_df = pd.DataFrame(new_cols)
    new_df = pd.concat([new_df, new_cols_df], axis=1)
    
    # sanity check
    for n_struct,old_list in roll_up_records.items():
        sum_old = df[old_list].sum(axis=1)
        sum_new = new_df[n_struct]
        assert sum(sum_old==sum_new) == len(df)
    
    
    
    # print(outfile)
    # print()
    assert os.path.abspath(outfile) != os.path.abspath(infile)
    new_df.to_csv(outfile)


def normalize_projection_columns_per_cell(input_df, projection_column_identifiers=['ipsi', 'contra']):
    """
    :param input_df:  input projection df
    :param projection_column_identifiers: list of identifiers for projection columns. i.e. strings that identify projection columns from metadata columns
    :return: normalized projection matrix
    """
    proj_cols = [c for c in input_df.columns if any([ider in c for ider in projection_column_identifiers])]
    input_df[proj_cols] = input_df[proj_cols].fillna(0)

    res = input_df[proj_cols].T / input_df[proj_cols].sum(axis=1)
    input_df[proj_cols] = res.T

    return input_df
