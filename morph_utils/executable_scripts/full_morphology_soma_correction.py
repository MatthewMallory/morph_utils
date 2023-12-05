import os
import pandas as pd
import argschema as ags
from multiprocessing import Pool
import matplotlib.pyplot as plt
from neuron_morphology.morphology import Morphology
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from morph_utils.graph_traversal import dfs_labeling, bfs_tree
from morph_utils.modifications import remove_duplicate_soma, re_root_morphology, sort_morph_ids, assign_soma_by_node_degree
from morph_utils.visuals import basic_morph_plot


class IO_Schema(ags.ArgSchema):
    swc_input_directory = ags.fields.InputDir(description='directory with micron resolution ccf registered files')
    swc_output_directory = ags.fields.OutputDir(description='directory with micron resolution ccf registered files')
    qc_image_output_directory = ags.fields.OutputDir(description="output projection csv")
    use_multiprocessing = ags.fields.Bool(default=False, description="weather to use cpu multiprocessing or not")
    allow_risk = ags.fields.Bool(default=False, description="allow algorithm to proceed and risk missing true soma")
    assign_soma_when_missing = ags.fields.Bool(default=True, description="if no soma, assign to most complex (# children)")


def correct_soma_in_swc_file(input_swc_path, output_swc_path, qc_img_outpath, allow_risk, assign_soma_when_missing):
    """
    This will fix the following:

    -- Incorrect soma was chosen by morphology.get_soma() because there are mutliple soma nodes in the swc file
    -- Children of the true soma are also labeled soma

    by assuming the node with the most children is the actual soma node. Children of soma nodes that are also labeled
    as type 1 (soma) will be traversed via depth first until a non soma node is reached. Then the traversed path will
    be retyped to whichever non soma node type was discovered.

    THIS CAN FAIL IF:
    Furcation nodes are represented by one root node for each furcation. For example, a bifurcation is represented as
    two uniquely rooted segments. Each root is a node located at the bifurcation coordinate whose parent is -1. This is
    a bug that shows up from editing swc files in certain softwares. If this is present, you should run the ivscc-sort
    plugin in vaa3d.

    BECAUSE:
    The reason this script will fail under those circumstances, is because the true soma node may be represented as
    individual uniquely rooted segments. For example if a neuron has 7  basal dendrite stems, but because of this
    bug those stems are represented as 7 uniquely rooted/individual segments, there is no way of identifying the
    true soma

    SUGGESTED:
    Run ivscc-sort in vaa3d This plugin can be run interatcively, or through command line. If through command line
    (batch processing) you will need soma nodes for each file to have the same node id. If this is not true for your
    files, you can run sort_morph_idsologies.py from  morph_utils/executable_scripts prior to vaa3d sort. This ensures the
    soma will be node id 1. After ivscc-sort in vaa3d you can then run this script.

    SUMMARY:
    This will return a dictionary of {swc_file:number root nodes} for a given input file if there are multiple root nodes
    This is considered risky because the true soma may be masked by the issue/bug described above. These files will
    be written to a csv called FailedFiles.csv in the qc image directory.

    :param input_swc_path: path to input swc file to be soma fixed
    :param output_swc_path: path to output swc file to be soma fixed
    :param qc_img_outpath: path to output qc before/after image
    :return: input_swc_path:number of root nodes or None if run to completion
    """
    morph = morphology_from_swc(input_swc_path)
    if morph.get_soma() is None:

        if assign_soma_when_missing:
            morph = assign_soma_by_node_degree(morph)
            if morph.get_soma() is None:
                return {"failed_filename": os.path.abspath(input_swc_path),
                        "fail_reason": "no soma node found by using number of children as guide"}

        else:
            return {"failed_filename": os.path.abspath(input_swc_path),
                    "fail_reason": "no soma node and assign_soma_when_missing assigned False"}

    root_nodes = [n for n in morph.nodes() if n['parent'] == -1]
    num_roots = len(root_nodes)
    retun_value = None
    if num_roots > 1:
        print("Too many root nodes ({}) in {}".format(num_roots, input_swc_path))
        if not allow_risk:
            return {"failed_filename": os.path.abspath(input_swc_path), "fail_reason": "too_many_roots ({}) while "
                                                                                       "allow_risk "
                                                                                       "flag=False".format(
                num_roots)}
        else:
            print("There are multiple root nodes in this file {}. ".format(input_swc_path))
            print("Because you designated --allow_risk True, you are trusting morph.get_soma()")
            print("The soma node indicated by morph.get_soma() = {}".format(morph.get_soma()))
            print("This file should be manually checked")
            retun_value = {"manuall_check_filename": os.path.abspath(input_swc_path), "number_root_nodes": num_roots}

    original_morph = morph.clone()
    original_soma = original_morph.get_soma()
    original_soma_id = original_soma['id']

    # Get num of duplicate soma nodes before fix
    soma_tuple = (original_soma['x'], original_soma['y'], original_soma['z'])
    num_nodes_at_soma_coord_before = [n for n in morph.nodes() if (n['x'], n['y'], n['z']) == soma_tuple]
    num_nodes_at_soma_coord_before = len(num_nodes_at_soma_coord_before) - 1  # -1 to account for the actual soma

    morph = remove_duplicate_soma(morph)
    soma_node = morph.get_soma()
    soma_tuple = (soma_node['x'], soma_node['y'], soma_node['z'])

    # Get num of duplicate soma nodes after fix
    num_nodes_at_soma_coord_after = [n for n in morph.nodes() if (n['x'], n['y'], n['z']) == soma_tuple]
    num_nodes_at_soma_coord_after = len(num_nodes_at_soma_coord_after) - 1  # -1 to account for the actual soma

    # Choose the most likely soma node
    soma_node_candidates = [n for n in morph.nodes() if n['type'] == 1]
    num_children_per_candidate = {n['id']: len(morph.get_children(n)) for n in soma_node_candidates}  # morph.nodes()}
    suspected_soma_id = \
    [k for k, v in num_children_per_candidate.items() if v == max(list(num_children_per_candidate.values()))][0]
    suspected_soma = morph.node_by_id(suspected_soma_id)

    if soma_node != suspected_soma:
        print("Found a replacement soma node")
        # for the segment that leads uptree from the new soma, we need to fix these so that they
        # are actually going down.
        morph = re_root_morphology(suspected_soma, morph)
        morph.node_by_id(suspected_soma_id)['parent'] = -1
        morph.node_by_id(suspected_soma_id)['type'] = 1
        morph = Morphology([n for n in morph.nodes()],
                           parent_id_cb=lambda x: x['parent'],
                           node_id_cb=lambda x: x['id'])
        soma_node = morph.node_by_id(suspected_soma_id)

    # Now we need to pull from outproblem one solution. neuron_morphology wont recognize our chosen
    # soma until after we relabel all other soma nodes to their appropriate label
    type_ct = {}
    soma_children = morph.get_children(soma_node)
    for comp in [1, 2, 3, 4]:
        type_ct[comp] = len([f for f in soma_children if f['type'] == comp])

    if type_ct[1] != 0:
        # There are children of the soma typed 1.
        # get the label of the somas grandchild and label the somas child as such
        for ch_no in soma_children:

            if ch_no['type'] == 1:
                grand_children = morph.get_children(ch_no)
                seg_down, _ = bfs_tree(ch_no, morph)

                # keep track of which nodes need to be re-typed (it may be more than just the soma's child)
                # as soon as we get to a non soma type we will stop.
                # This assumes a segment is not all labeled as type = 1
                nodes_to_retype = []
                for no in seg_down:
                    if no['type'] != 1:
                        new_label = no['type']
                        break
                    else:
                        nodes_to_retype.append(no)

                for node in nodes_to_retype:
                    morph.node_by_id(node['id'])['type'] = new_label

        # make a new morph object
        morph = Morphology(morph.nodes(),
                           parent_id_cb=lambda x: x['parent'],
                           node_id_cb=lambda x: x['id'])
        # re calculate the type count dict
        new_type_ct = {}
        soma_children = morph.get_children(soma_node)
        for comp in [1, 2, 3, 4]:
            new_type_ct[comp] = len([f for f in soma_children if f['type'] == comp])

        # make sure  changes worked
        assert new_type_ct[1] == 0

    morph = sort_morph_ids(morph, soma_node=soma_node)
    final_soma_id = morph.get_soma()['id']
    morphology_to_swc(morph, output_swc_path)

    fig, ax = plt.subplots(1, 2)
    basic_morph_plot(original_morph, ax[0], title="Original. SomaID {}".format(original_soma_id))
    basic_morph_plot(morph, ax[1], title="Fixed. SomaID {}".format(final_soma_id))

    for i in [0, 1]:
        buff = 20
        ax[i].set_xlim(soma_node['x'] - buff, soma_node['x'] + buff)
        ax[i].set_ylim(soma_node['y'] - buff, soma_node['y'] + buff)
        ax[i].set_aspect('equal')
        ax[i].legend()
    fig.set_size_inches(14, 6)

    fig.savefig(qc_img_outpath, dpi=300, bbox_inches='tight')
    plt.clf()
    return retun_value



def main(swc_input_directory, swc_output_directory, qc_image_output_directory, use_multiprocessing, allow_risk,
         assign_soma_when_missing,**kwargs):
    swc_input_directory = os.path.abspath(swc_input_directory)
    swc_output_directory = os.path.abspath(swc_output_directory)
    qc_image_output_directory = os.path.abspath(qc_image_output_directory)

    parallel_inputs = []
    reslist = []
    for swc_fn in [f for f in os.listdir(swc_input_directory) if f.endswith(".swc")]:
        swc_src_path = os.path.join(swc_input_directory, swc_fn)
        swc_dest_path = os.path.join(swc_output_directory, swc_fn)
        qc_img_path = os.path.join(qc_image_output_directory, swc_fn.replace(".swc", ".png"))


        if use_multiprocessing:
            parallel_inputs.append((swc_src_path, swc_dest_path, qc_img_path,allow_risk,assign_soma_when_missing))
        else:
            reslist.append(correct_soma_in_swc_file(swc_src_path, swc_dest_path, qc_img_path,allow_risk,assign_soma_when_missing))
    if use_multiprocessing:
        p = Pool()
        reslist = p.starmap(correct_soma_in_swc_file, parallel_inputs)

    reslist = [f for f in reslist if f != None]
    if reslist != []:
        df_ofile = os.path.join(qc_image_output_directory, "Failed_Files.csv")
        resdf = pd.DataFrame.from_records(reslist)
        resdf.to_csv(df_ofile)

    print(len(reslist))

def console_script():
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)

if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(**module.args)
