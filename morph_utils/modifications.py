import os
import numpy as np
from collections import deque
from scipy.spatial.distance import euclidean
from neuron_morphology.morphology import Morphology
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from neuron_morphology.transforms.affine_transform import AffineTransform
from morph_utils.graph_traversal import dfs_labeling, bfs_tree, get_path_to_root
from morph_utils.query import query_for_z_resolution
from morph_utils.measurements import dist_bwn_nodes
from scipy import interpolate
from copy import copy


def prune_tree(morphology,num_node_thresh, node_types=[1,2,3,4] ):
    """will prune any segments in the tree that are shorter than a given length threhsold, 
    where length is measured in number of nodes

    Args:
        morphology (NeuronMorphology.morphology): input morphology
        num_node_thresh (int): pruning threshold
        node_types (list): node types to consider when pruning 

    Returns:
        NeuronMorphology.morphology: pruned morphology
    """
    nodes_to_remove = set()
    prune_count = 0
    bifur_nodes = [n for n in morphology.nodes() if (len(morphology.get_children(n)) > 1) and (n['type'] in node_types) ]
    for bif_node in bifur_nodes:
        children = morphology.get_children(bif_node)
        for child in children:
            child_remove_nodes, child_seg_length = bfs_tree(child, morphology)
            if child_seg_length < num_node_thresh:
                prune_count += 1
                [nodes_to_remove.add(n['id']) for n in child_remove_nodes]

    soma = morphology.get_soma()
    soma_children = [ch for ch in morphology.get_children(soma) if ch['type'] in node_types]
    soma_children_ids = [n['id'] for n in soma_children]
    root_nodes = soma_children + [n for n in morphology.nodes() if n['parent']==-1 and (n['id'] not in soma_children_ids) and (n['id']!=soma['id']) ]
    root_nodes = [r for r in root_nodes if r['type'] in node_types]
    for root in root_nodes:
        down_tree,down_tree_n = bfs_tree(root,morphology)
        if down_tree_n < num_node_thresh:
            prune_count += 1
            [nodes_to_remove.add(n['id']) for n in down_tree]
            
    keeping_nodes = [n for n in morphology.nodes() if n['id'] not in nodes_to_remove]
    pruned_morph = Morphology(
        keeping_nodes,
        node_id_cb=lambda node: node['id'],
        parent_id_cb=lambda node: node['parent'])

    return pruned_morph

def resample_3d_points(points, spacing):
    """
    Resample points at a given spacing. Will include the first and last points provided in points. 
    

    Args:
        points (nxm np array): must have at least two data points. 
        spacing (float): desired spacing to resample data points. Must be positive number greater than zero

    Returns:
        np array: resampled numpy array
    """
    # Extract x, y, and z coordinates from the points array
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    # Calculate the cumulative distance along the curve
    cumulative_distance = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2))
    cumulative_distance = np.insert(cumulative_distance, 0, 0)  # Add a zero at the beginning

    # Create an interpolation function for each coordinate
    interpolate_x = interpolate.interp1d(cumulative_distance, x, kind='linear', fill_value="extrapolate")
    interpolate_y = interpolate.interp1d(cumulative_distance, y, kind='linear', fill_value="extrapolate")
    interpolate_z = interpolate.interp1d(cumulative_distance, z, kind='linear', fill_value="extrapolate")

    # Create a new set of distances with the desired spacing
    new_distances = np.arange(0, cumulative_distance[-1], spacing)

    # Interpolate the coordinates at the new distances
    new_x = interpolate_x(new_distances)
    new_y = interpolate_y(new_distances)
    new_z = interpolate_z(new_distances)

    # Create the resampled array
    resampled_array = np.column_stack((new_x, new_y, new_z))

    if len(resampled_array)==1:
        # the desired spacing was larger than the max extend between the start and finish nodes of this segment,
        # so just return those
        return np.vstack([points[0],points[-1]])
    
    return np.vstack([resampled_array, points[-1]])


def resample_morphology(morph, spacing_size):
    """
    Will resample the spacing between ancestor-descendant irreducible node pairs. In this function, 
    we consider irreducible nodes to be leaf nodes, branch nodes and the immediate children of the soma.
    The immediate children of the soma are considered irreducible so that resampling does not 
    occur between the soma and it's first descendant because this space should be occupied by 
    the somas' radius. 
    
    If the spacing_size is larger than the distance between a given pair of ancestor-descendant 
    irreducible nodes, only the irreducible nodes will remain in the morphology. 
    
    Args:
        morph (neuron_morphology.Morphology): input morphology
        spacing_size (float): desired spacing between nodes. 
    """

    # iterate over roots so this can handle autotrace cells that have multiple roots (disconnected segments)
    roots = [n for n in morph.nodes() if n['parent']==-1 and n['type']==1]
    roots = roots + [n for n in morph.nodes() if n['parent']==-1 and n['type']!=1]
    new_nodes = []
    node_ct = 1

    old_irr_id_to_new_irr_id_dict = {}
    for root in roots:
        new_root = copy(root)
        new_root['id'] = node_ct
        new_nodes.append(new_root)
        old_irr_id_to_new_irr_id_dict[root['id']]=new_root['id']
        this_roots_children = morph.get_children(root)
        node_ct+=1
          
        for child in this_roots_children:
                
            new_child = copy(child)
            new_child['parent'] = new_root['id']
            new_child['id'] = node_ct
            new_nodes.append(new_child)
            old_irr_id_to_new_irr_id_dict[child['id']]=new_child['id']
                    
            # get a list of lists where each sublist is a list of nodes connecting two irreducible 
            # nodes. The first node will be the upstream irreducible node and the last node will be the
            # descendant irreducible node.
            this_list = []
            irreducible_segments = []
            queue = deque([child])
            seen_ids = set()
            while len(queue) > 0:
                
                current_node = queue.popleft() 
                
                parent = morph.node_by_id(current_node['parent'])
                siblings = morph.get_children(parent)
                if len(siblings)>1 and parent['id'] in seen_ids:
                    if parent['id']!=morph.get_soma()['id']:
                        this_list.append(parent)
                    
                this_list.append(current_node)
                seen_ids.update([current_node['id']])
                children_list = morph.get_children(current_node)
                if len(children_list)!=1:
                    irreducible_segments.append(this_list)
                    this_list = []
                
                for ch_no in children_list:
                    queue.appendleft(ch_no)


            already_seen_irr_node_ids = set()
            for sublist in irreducible_segments:
                # get the ancestor and descendant irreducible nodes
                irr_node_1 = sublist[0]
                irr_node_2 = sublist[-1]
                
                dist_between_irr_nodes = dist_bwn_nodes(irr_node_1, irr_node_2)
                segment_arr = np.array([[n['x'],n['y'],n['z']] for n in sublist])
                if (segment_arr.shape[0] > 2) or (dist_between_irr_nodes>spacing_size) :
                    # resample the segment when there are multiple nodes, or the
                    # space between the nodes is greater than sampling size
                    segment_arr_resamp = resample_3d_points(segment_arr, spacing_size)
                    reducible_arr = segment_arr_resamp[1:-1]
                else:
                    reducible_arr = []
                    
                # determine what node id the first reducible node should point to                
                if irr_node_1['id'] in already_seen_irr_node_ids:
                    
                    # we have recursed back up the tree to go down a different branch
                    # we need our first reducible node to point to the new id assigned 
                    # to irr_node_1 when we saw it first in a previous iteration
                    red_1_parent_id = old_irr_id_to_new_irr_id_dict[irr_node_1['id']]
                    
                else:
                    # we have not recured, still moving down a segment in dfs.
                    red_1_parent_id = node_ct
                
                if len(reducible_arr) != 0:
                        
                    # setup an equivalent index and step size to make sure
                    # that we are using approriate radius information 
                    # in the resampled morphology. 
                    step_size = len(segment_arr)/len(segment_arr_resamp)
                    equiv_idx = 0
                    for new_coord_ct, c in enumerate(reducible_arr):

                        equiv_node = sublist[equiv_idx]
                        node_ct+=1
                        if new_coord_ct==0:
                            parent_id = red_1_parent_id
                        else:
                            parent_id = node_ct-1
                        new_node = {
                            "x":c[0],
                            "y":c[1],
                            "z":c[2],
                            'id':node_ct,
                            "type":equiv_node['type'],
                            "radius":equiv_node["radius"],
                            "parent":parent_id,
                        }
                        new_nodes.append(new_node)
                        
                        equiv_idx+=step_size
                        equiv_idx = np.math.floor(equiv_idx)
                                
                            
                node_ct+=1
                new_node_2 = copy(irr_node_2)
                if len(reducible_arr)==0:
                    new_node_2['parent'] = red_1_parent_id
                else:
                    new_node_2['parent'] = node_ct-1
                new_node_2['id']=node_ct
                new_nodes.append(new_node_2)
                
                old_irr_id_to_new_irr_id_dict[irr_node_2['id']] = new_node_2['id']
                
                already_seen_irr_node_ids.update([irr_node_1['id']])

            node_ct+=1    

        resampled_morph = Morphology(new_nodes,
                parent_id_cb=lambda x:x['parent'],
                node_id_cb=lambda x:x['id'])
        
    return resampled_morph


def generate_irreducible_morph(morph):
    """
    Will generate an irreducible morphology object. The only remaining nodes will be roots,
    branches and tip nodes.

    :param morph: neuron_morphology Morphology object
    :return: neuron_morphology Morphology object
    """
    morph = morph.clone()

    irreducible_nodes = [n for n in morph.nodes() if
                         (len(morph.get_children(n)) > 1) or (len(morph.get_children(n)) == 0) or (n['parent'] == -1)]
    soma = morph.get_soma()
    if not soma:
        soma_list = [n for n in morph.nodes() if n['parent'] == -1]
        if len(soma_list) != 1:
            print("Invalid Number of somas (0 or >1)")
            return None
        else:
            soma = soma_list[0]
    if soma not in irreducible_nodes and soma:
        irreducible_nodes.append(morph.get_soma())

    leaves = [n for n in morph.nodes() if len(morph.get_children(n)) == 0]
    irreducible_nodes_with_topology = []
    # need to re-assign parent child relationship for only irreducible nodes
    for leaf_no in leaves:
        path_to_root = get_path_to_root(leaf_no, morph)

        if leaf_no not in path_to_root:
            path_to_root.insert(0, leaf_no)

        irreducible_nodes_in_path = [n for n in path_to_root if n in irreducible_nodes]

        for i in range(0, len(irreducible_nodes_in_path) - 1):
            this_no = irreducible_nodes_in_path[i]
            next_node_up = irreducible_nodes_in_path[i + 1]

            this_no['parent'] = next_node_up['id']
            irreducible_nodes_with_topology.append(this_no)

        # add root 
        next_node_up['type'] = 1
        irreducible_nodes_with_topology.append(next_node_up)

    morph_irreducible = Morphology(irreducible_nodes_with_topology,
                                   parent_id_cb=lambda x: x['parent'],
                                   node_id_cb=lambda x: x['id'])

    return morph_irreducible


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


def normalize_position(morph): 
    """
    Shift morphology position so the soma is at 0,0,0

    :param morph: neuron_morphology Morphology object
    :return: noramlized neuron_morphology Morphology object 
    """
    soma = morph.get_soma()
    trans_list = [1,0,0, 0,1,0, 0,0,1, -soma['x'],-soma['y'],-soma['z']]
    translate_transform= AffineTransform.from_list(trans_list)
    moved_morph = translate_transform.transform_morphology(morph) # if you need the original object to remain unchanged do morph.clone()

    return moved_morph

def convert_pixel_to_um(morph, specimen_id):
    """
    Convert morphology units from pixel to micron. 

    :param morph: neuron_morphology Morphology object in pixel units
    :param specimen_id: cell specimen id
    :return: neuron_morphology Morphology object in micron units
    """
    anisotropy_value = query_for_z_resolution(specimen_id)
    scale_list = [0.1144,0,0,  0,0.1144,0,   0,0,anisotropy_value,   0,0,0] 
    scale_transform = AffineTransform.from_list(scale_list)
    scaled_morph = scale_transform.transform_morphology(morph) # if you need the original object to remain unchanged do morph.clone()

    return scaled_morph