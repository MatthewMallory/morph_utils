import os
import argschema as ags
from multiprocessing import Pool
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from morph_utils.modifications import remove_duplicate_soma, sort_morph_ids


class IO_Schema(ags.ArgSchema):
    swc_input_directory = ags.fields.InputDir(description='directory with micron resolution ccf registered files')
    swc_output_directory = ags.fields.OutputDir(description='directory with micron resolution ccf registered files')
    use_multiprocessing = ags.fields.Bool(default=False, description="weather to use cpu multiprocessing or not")


def sort_swc_file(input_swc_path, output_swc_path):
    """
    This will run duplicate (stacked) soma node remover and then sort_morph_ids function across directory of swc file.
    sort_morph_ids will not assign soma id to 1 if there are multiple somas stacked on top of eachother so that is why we
    run that first.s
    :param input_swc_path: path to input swc file to be soma fixed
    :param output_swc_path: path to output swc file to be soma fixed
    :return: None/input_swc_path if failed
    """

    morph = morphology_from_swc(input_swc_path)
    morph = remove_duplicate_soma(morph)
    morph = sort_morph_ids(morph)
    morphology_to_swc(morph, output_swc_path)


def main(swc_input_directory, swc_output_directory, use_multiprocessing, **kwargs):
    parallel_inputs = []
    reslist = []
    for swc_fn in os.listdir(swc_input_directory):
        swc_src_path = os.path.join(swc_input_directory, swc_fn)
        swc_dest_path = os.path.join(swc_output_directory, swc_fn)

        if use_multiprocessing:
            parallel_inputs.append((swc_src_path, swc_dest_path))
        else:
            reslist.append(sort_swc_file(swc_src_path, swc_dest_path))
    if use_multiprocessing:
        p = Pool()
        reslist = p.starmap(sort_swc_file, parallel_inputs)

    # reslist = [f for f in reslist if f != None]
    # if reslist != []:
    #     df_ofile = os.path.join(swc_output_directory, "Failed_Files.csv")
    #     resdf = pd.DataFrame.from_records(reslist)
    #     resdf.to_csv(df_ofile)
    #
    # print(len(reslist))

def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(module.args)


if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
