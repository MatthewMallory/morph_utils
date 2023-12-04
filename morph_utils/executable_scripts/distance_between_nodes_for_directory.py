import os
import pandas as pd
import argschema as ags
import itertools
from morph_utils.measurements import node_distance_between_morphs
from multiprocessing import Pool
import numpy as np


class IO_Schema(ags.ArgSchema):
    swc_input_directory = ags.fields.InputDir(description='directory with swc files')
    output_file = ags.fields.OutputFile(descripion='output csv with distances between files')
    compartment_types = ags.fields.List(default=[3, 4], cls_or_instance=ags.fields.Int)
    use_multiprocessing = ags.fields.Boolean(default=True)


def main(swc_input_directory, output_file, compartment_types, use_multiprocessing, **kwargs):
    """
    Will create csv reporting mean distance between nodes of a certain type (axon/basal/apical dendrite)
    for all combinations of swc files in directory. This will leverage parallel processing to speed things up.

    This is particularly useful in trying to identify duplicate/very similar swc files, but may also be used as a
    rudimentary asymmetric similarity metric for comparing neurons.
    """
    parallel_inputs = []
    swc_input_directory = os.path.abspath(swc_input_directory)
    file_list = [os.path.join(swc_input_directory,f) for f in os.listdir(swc_input_directory)]
    all_combinations = list(itertools.combinations(file_list, 2))
    # all_combinations = [c for c in all_combinations] # if c[0] != c[1]]

    print("{} Comparisons to analyze".format(len(all_combinations)))
    reslist = []
    for combo in all_combinations:
        file_1 = os.path.join(swc_input_directory, combo[0])
        file_2 = os.path.join(swc_input_directory, combo[1])

        if use_multiprocessing:
            parallel_inputs.append((file_1, file_2, compartment_types))
        else:
            reslist.append(node_distance_between_morphs(file_1, file_2, compartment_types))

    if use_multiprocessing:
        p = Pool()
        reslist = p.starmap(node_distance_between_morphs, parallel_inputs)

    # reslist is a list of dictionaries
    resulting_df = pd.DataFrame.from_records(reslist)
    resulting_df.to_csv(output_file.replace(".csv", "_long.csv"))

    # Build distance matrix from this long format
    distance_matrix = np.zeros((len(file_list), len(file_list)))
    for combo in all_combinations:
        file_1 = os.path.join(swc_input_directory, combo[0])
        file_2 = os.path.join(swc_input_directory, combo[1])

        idx_1, idx_2 = file_list.index(file_1), file_list.index(file_2)
        if file_1 == file_2:
            distance_matrix[idx_1, idx_2] = 1
        else:
            this_row = resulting_df[(resulting_df['file_1'] == file_1) & (resulting_df['file_2'] == file_2)]
            forward_score = this_row['forward_distance'].iloc[0]
            reverse_score = this_row['reverse_distance'].iloc[0]

            distance_matrix[idx_1, idx_2] = forward_score
            distance_matrix[idx_2, idx_1] = reverse_score

    distance_df = pd.DataFrame(distance_matrix, columns=file_list, index=file_list)
    distance_df.to_csv(output_file)


def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)


if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
