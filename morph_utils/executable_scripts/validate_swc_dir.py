import os
import ast
import argschema as ags
from multiprocessing import Pool
import pandas as pd
from morph_utils.validation import ivscc_validate_morph


class IO_Schema(ags.ArgSchema):
    swc_input_directory = ags.fields.InputDir(description='directory with swc files to validate')
    report_csv = ags.fields.OutputFile(description='directory with micron resolution ccf registered files')
    soma_child_distance_threshold = ags.fields.Int(default=50,
                                                   descreiption='max distance a somas child may be from soma')
    use_multiprocessing = ags.fields.Bool(default=True, description="weather to use cpu multiprocessing or not")


def main(swc_input_directory, soma_child_distance_threshold, report_csv, use_multiprocessing,
         **kwargs):
    print("Soma-Child Distance Threshold: {}".format(soma_child_distance_threshold))
    print("Use Multiprocessing: {}".format(use_multiprocessing))
    parallel_inputs = []
    reslist = []
    for swc_fn in os.listdir(swc_input_directory):
        swc_src_path = os.path.join(swc_input_directory, swc_fn)
        if use_multiprocessing:
            parallel_inputs.append((swc_src_path, soma_child_distance_threshold))
        else:
            reslist.append(ivscc_validate_morph(swc_src_path, soma_child_distance_threshold))
    if use_multiprocessing:
        p = Pool()
        reslist = p.starmap(ivscc_validate_morph, parallel_inputs)

    res_df = pd.DataFrame.from_records(reslist)
    # res_df['error_list'] = res_df['error_list'].apply(lambda x: ast.literal_eval(x))
    max_num_errs = res_df['error_list'].apply(lambda x: len(x)).max()
    err_cols = ["Error_{}".format(i) for i in range(max_num_errs)]

    records = []
    for idx, row in res_df.iterrows():

        sp_dict = {"full_swc_path": row.file_name,
                   "swc_file_name": os.path.basename(row.file_name)}
        for ec in err_cols:
            sp_dict[ec] = None

        err_lst = row.error_list
        for ct, sp_error in enumerate(err_lst):
            sp_dict['Error_{}'.format(ct)] = sp_error

        if len(err_lst) != 0:
            records.append(sp_dict)

    final_qc_df = pd.DataFrame.from_records(records)

    final_qc_df.to_csv(report_csv)


def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)


if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
