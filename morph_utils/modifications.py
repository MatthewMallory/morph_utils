import os
import numpy as np
from collections import deque
from scipy.spatial.distance import euclidean
from neuron_morphology.morphology import Morphology
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from morph_utils.graph_traversal import dfs_labeling, bfs_tree, get_path_to_root


def assign_soma_by_node_degree(morphology, num_children_threshold=2):
    """
    Will assign soma to the node that has the most children. This will NOT remove duplicate soma nodes,
    only assign the node with highest degree as soma.

    :param num_children_threshold: the minimum number of children a true soma node will have
    :param morphology: neuron_morphology Morphology object
    :return: neuron_morphology Morphology object
    """

    soma_types = [n for n in morphology.nodes() if n['type'] == 1]
    if len(soma_types) != 1:

        num_children_per_node = {n['id']: len(morphology.get_children(n)) for n in morphology.nodes()}
        max_num_children = max(list(num_children_per_node.values()))
        if max_num_children >= num_children_threshold:
            no_ids = [k for k, v in num_children_per_node.items() if v == max_num_children]

            if len(no_ids) > 1:
                # find which node is closest to the morphology centroid
                coords = np.array([[n['x'], n['y'], n['z']] for n in morphology.nodes()])
                center = np.mean(coords, axis=0)
                centroid = (center[0], center[1], center[2])

                min_dist_to_centroid = np.inf
                chosen_node = no_ids[0]
                for no_id in no_ids:
                    this_no = morphology.node_by_id(no_id)
                    this_no_coord = (this_no['x'], this_no['y'], this_no['z'])
                    this_dist_to_centroid = euclidean(this_no_coord, centroid)

                    if this_dist_to_centroid < min_dist_to_centroid:
                        min_dist_to_centroid = this_dist_to_centroid
                        chosen_node = no_id
            else:
                chosen_node = no_ids[0]
                print(
                    "Choosing new soma based on num children. There are {} nodes with max value of {} children".format(
                        len(no_ids),
                        max_num_children))

            morphology.node_by_id(chosen_node)['type'] = 1
            keeping_nodes = morphology.nodes()

            new_morph = Morphology(keeping_nodes,
                                   node_id_cb=lambda x: x['id'],
                                   parent_id_cb=lambda x: x['parent'])
            print("New Morphs Soma = {}".format(new_morph.get_soma()))

            return new_morph
        else:
            print("There are no nodes in the morphology that have at least {} children".format(num_children_threshold))
            return morphology

    else:
        return morphology


def remove_duplicate_soma(morphology, soma=None):
    """
    Will remove nodes that are at the same coordinate as the soma node (regardless of node type) where
    soma node defaults to that found by morphology.get_soma(). Any child of  a duplicate soma node
    will be adopted by the soma node.

    This will NOT remove nodes in the morphology that are type 1, but at a different coordinate than the soma

    :param morphology: neuron_morphology Morphology object
    :return: neuron_morphology Morphology object
    """
    morphology = morphology.clone()

    if soma is None:
        soma = morphology.get_soma()

    duplicate_somas = [n for n in morphology.nodes() if (n['x'], n['y'], n['z']) == (soma['x'], soma['y'], soma['z'])]
    duplicate_somas.remove(soma)

    if duplicate_somas == []:
        # no duplicate somas, just check that soma's parent is -1 and it's id is 1
        if soma is not None:
            if soma['id'] != 1:
                morphology = sort_morph_ids(morphology, soma_node=soma)
                soma = morphology.get_soma()

            if soma['parent'] != -1:
                print("UH SOMA'S PARENT IS NOT = -1")
        else:
            print("No Soma?")
        return morphology

    duplicate_soma_ids = [n['id'] for n in duplicate_somas]

    # find actual children of soma. These are the children of our duplicate soma nodes that aren't also in duplicate soma
    # node, but also
    children_of_soma = [n for n in morphology.nodes() if
                        (n['parent'] in duplicate_soma_ids) and (n not in duplicate_somas)]

    # assign their parent to the chosen soma node
    for no in children_of_soma:
        morphology.node_by_id(no['id'])['parent'] = soma['id']

    # make sure the somas parent is -1 and type is 1
    morphology.node_by_id(soma['id'])['parent'] = -1
    morphology.node_by_id(soma['id'])['type'] = 1

    # create new morphology
    keeping_nodes = [n for n in morphology.nodes() if n['id'] not in duplicate_soma_ids]
    new_morph = Morphology(keeping_nodes,
                           node_id_cb=lambda x: x['id'],
                           parent_id_cb=lambda x: x['parent'])
    new_soma = new_morph.get_soma()

    # sort if needed so that soma id is = 1
    new_soma = new_morph.get_soma()
    if new_soma['id'] != 1:
        new_morph = sort_morph_ids(new_morph)

    return new_morph


def sort_morph_ids(morph, soma_node=None, specimen_id=None, **kwargs):
    """
    Will sort a moprhology so that node id ascends from soma in a depth  first order. Will assure that the
    soma id is equal to 1

    TODO update so that we are not using IO operations and just creating a new morphology...

    :param morph: neuron_morphology Morphology object
    :param soma_node: a soma node (dictionary) from neuron_morphology Morphology object
    :param specimen_id: not required, used for temporary file naming
    :return:
    """
    if specimen_id is None:
        specimen_id = np.random.randint(0, 100000000)

    unsorted_swc_path = '{}_temp_sorting.swc'.format(specimen_id)
    sorted_swc_path = '{}_temp_sorted.swc'.format(specimen_id)
    morphology_to_swc(morph, unsorted_swc_path)
    unordered_swc_info = {}
    with open(unsorted_swc_path, 'r') as f:
        for l in f:
            if '#' not in l:
                no_id = int(l.split(' ')[0])
                parent_id = l.split()[-1]
                children_list = morph.get_children(morph.node_by_id(no_id))
                unordered_swc_info[no_id] = l

    new_node_ids = {}
    start_label = 1
    if soma_node is None:
        soma_node = morph.get_soma()
    #         root_node_list = morph.get_roots()
    #     else:
    #         root_node_list = [n for n in morph.nodes() if n['parent']==-1 and n['type']==1]#morph.get_roots()

    root_node_list = morph.get_roots() + [soma_node]
    unique_root_ids = set([n['id'] for n in root_node_list])
    root_node_list = [morph.node_by_id(i) for i in unique_root_ids]

    root_node_list.remove(soma_node)

    # Start with soma so its node id is one
    seg_len = dfs_labeling(soma_node, start_label, new_node_ids, morph)
    start_label += seg_len

    for root in root_node_list:
        seg_len = dfs_labeling(root, start_label, new_node_ids, morph)
        start_label += seg_len

    new_output_dict = {}
    # with open(sorted_swc_path,"w") as f2:
    for old_id, old_line in unordered_swc_info.items():
        new_id = new_node_ids[old_id]
        old_parent = int(old_line.split()[-1])
        if old_parent == -1:
            new_parent = -1
        else:
            new_parent = new_node_ids[old_parent]

        new_line_list = [str(new_id)] + old_line.split(' ')[1:-1] + ['{}\n'.format(new_parent)]
        new_line = " ".join(new_line_list)
        new_output_dict[new_id] = new_line
        # f2.write(new_line)

    with open(sorted_swc_path, "w") as f2:
        for k in sorted(list(new_output_dict.keys())):
            new_write_line = new_output_dict[k]
            f2.write(new_write_line)

    sorted_morph = morphology_from_swc(sorted_swc_path)
    os.remove(sorted_swc_path)
    os.remove(unsorted_swc_path)
    return sorted_morph


def re_structure_segment(morphology, new_root_node, new_roots_parent=-1, overwrite_soma=True):
    """
    Don't think this is actively being used anywhere. Holding on to it for now until this package is more finalized

    :param morphology:
    :param new_root_node:
    :param new_roots_parent:
    :param overwrite_soma:
    :return:
    """
    morphology = morphology.clone()
    path_up = get_path_to_root(new_root_node, morphology)
    path_up = [n for n in path_up if n != new_root_node]
    current_root = path_up[-1]

    path_down, _ = bfs_tree(new_root_node, morphology)
    path_down = [n for n in path_down if n != new_root_node]

    if (not overwrite_soma) and (current_root['type'] == 1):
        print("You are trying to re-root a segment that has a node of type 1 as it's curent root")
        return morphology

    ct = -1
    for no_up in path_up:
        ct += 1
        # only take into consideration children that are in our direct path from new root to current root.
        children = [n for n in morphology.get_children(no_up) if n in path_up]

        if (ct == 0) and (children == []):
            # in this scenario we want to make sure our new root node is the first nodes parent. For whatever
            # vaa3d reason this first node in path_up has no children? To follow the logic of merging all soma
            # nodes later in this script, we need them to all be roots
            children = [new_root_node]

        assert len(children) == 1
        future_parent = children[0]
        morphology.node_by_id(no_up['id'])['parent'] = future_parent['id']

    morphology.node_by_id(new_root_node['id'])['parent'] = new_roots_parent

    new_nodes = [n for n in morphology.nodes()]
    new_morph = Morphology(new_nodes,
                           node_id_cb=lambda x: x['id'],
                           parent_id_cb=lambda x: x['parent'])
    return new_morph


def strip_compartment_from_morph(morph, compartment_to_strip):
    """
    remove all nodes of a certain type from a morphology

    :param morph: a neuron_morphology Morphology object
    :param compartment_to_strip: list of compartment types to remove [e.g. compartment = [3,4] would leave only soma and axon nodes]
    :return: neuron_morphology Morphology object
    """
    nodes = [n for n in morph.nodes() if n['type'] not in compartment_to_strip]
    axon_strip_morph = Morphology(nodes,
                                  parent_id_cb=lambda x: x['parent'],
                                  node_id_cb=lambda x: x['id'])

    return axon_strip_morph


def check_morph_for_segment_restructuring(morph):
    """
    This function will check the roots of all disconnected segments by visiting each root node and ensure
    that the closest leaf node to soma is the root. Particularly useful for autotrace processing

    :param morph: neuron_morphology Morphology object
    :return: morphology, Bool: Will return the unedited or edited morphology
    if needed, and True/False to represent re-structuring changes were made or not
    """
    morph = morph.clone()
    soma = morph.get_soma()
    if soma is not None:
        soma_coord = (soma['x'], soma['y'], soma['z'])
        roots = [n for n in morph.get_roots() if n != soma]  # + morph.get_children(soma)

        changes = False
        for root_node in roots:
            root_no_coord = (root_node['x'], root_node['y'], root_node['z'])
            root_dist_to_soma = euclidean(root_no_coord, soma_coord)
            seg_down, _ = bfs_tree(root_node, morph)
            leaf_nodes = [n for n in seg_down if morph.get_children(n) == []]

            closest_dist = root_dist_to_soma
            closest_node = root_node
            for leaf_no in leaf_nodes:
                leaf_no_coord = (leaf_no['x'], leaf_no['y'], leaf_no['z'])
                leaf_dist_to_soma = euclidean(leaf_no_coord, soma_coord)

                if leaf_dist_to_soma < closest_dist:
                    changes = True
                    closest_dist = leaf_dist_to_soma
                    closest_node = leaf_no

            if closest_node != root_node:
                morph = re_root_morphology(new_start_node=closest_node,
                                           morphology=morph)

        return morph, changes

    else:
        return morph, False


def re_root_morphology(new_start_node, morphology):
    """
    Will reorganize nodes so that the new_start_node becomes the root and the old root becomes a leaf node. Particularly
    useful in "flipping" parent child relationship direction in a disconnected auto-trace segment.

    formerly called restructure_disconnected_segment

    :param new_start_node: node that will now become root
    :param morphology: neuron_morphology Morphology object
    :return: re-rooted morphology
    """
    queue = deque()
    queue.append(new_start_node)

    new_parent_dict = {}
    visited_ids = [-1]
    while len(queue) > 0:
        this_node = queue.popleft()
        new_parent_dict[this_node['id']] = visited_ids[-1]
        visited_ids.append(this_node['id'])
        parent_id = this_node['parent']
        if parent_id != -1:
            parent_node = morphology.node_by_id(parent_id)
            queue.append(parent_node)

    for no_id, new_parent_id in new_parent_dict.items():
        morphology.node_by_id(no_id)['parent'] = new_parent_id

    new_morph = Morphology([n for n in morphology.nodes()],
                           node_id_cb=lambda x: x['id'],
                           parent_id_cb=lambda x: x['parent'])

    return new_morph
