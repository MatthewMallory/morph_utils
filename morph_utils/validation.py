from scipy.spatial.distance import euclidean
from collections import defaultdict
from morph_utils.graph_traversal import dfs_loop_check
from neuron_morphology.swc_io import morphology_from_swc
import os


def ivscc_validate_morph(input_swc, distance_threshold=50, expected_types=[1, 2, 3, 4]):
    """
    Standard protocol for validating an swc file as per AIBS IVSCC pipeline. This will:
    1. Make sure a soma node exists
    2. Make sure the soma is not > distance_threshold from it's children
    3. The immediate children of the soma are not furcation nodes
    4. There is only one root node (i.e. one node who has parent = -1)
    5. There is only 1 place of axon origination
    6. All nodes parents' exist in the tree (with the exception of a parent=-1, this will go towards the root node count)
    7. There is only one node of type = 1
    8. There are no nodes that share identical x,y,z coordinates
    9. All apical/basal dendrite nodes have parent of either soma or apical/basal dendrite respectively
    10. All axon nodes have parent of either axon, soma, or basal dendrite.
    11. There are no loops
    12. The soma ID is 1
    13. TODO Check if node IDs are sorted

    :param input_swc: path to swc file; str
    :param distance_threshold: maximum distance a valid child may be from soma
    :param expected_types: expected node types
    :return: error_list: list of all errors encountered; list
    """
    
    meta_records = {"errors":[]}
    node_list = []
    axon_origin_node_list = []
    root_node_list = [] 
    
    morph = morphology_from_swc(input_swc)

    nodes_to_qc = morph.nodes()
    all_node_ids = [n['id'] for n in nodes_to_qc]

    error_list = []
    axon_origins = 0
    number_of_roots = 0
    number_nodes_per_type = defaultdict(int)
    number_of_nodes_at_each_coord = defaultdict(int)
    node_coord_to_id = {}
    soma = morph.get_soma()
    if soma is None:
        error_list.append("No Soma Found")
        node_list.append({})
        error_dict = {'error':"No Soma Found", "Nodes":[]}
        meta_records['errors'].append(error_dict)
    else:
        soma_coord = (soma['x'], soma['y'], soma['z'])

    for no in nodes_to_qc:
        cur_node_type = no['type']
        parent_id = no['parent']
        cur_node_coord = (no['x'], no['y'], no['z'])

        number_nodes_per_type[cur_node_type] += 1
        number_of_nodes_at_each_coord[cur_node_coord] += 1
        node_coord_to_id[cur_node_coord] = no['id']
        
        if parent_id != -1:
            if parent_id in all_node_ids:

                if cur_node_type not in expected_types:
                    error_list.append("Unexpected Node Type Found: {}".format(cur_node_type))
                    node_list.append(no)
                    error_dict = {'error':"Unexpected Node Type Found: {}".format(cur_node_type), 
                                  "Nodes":[no]}
                    meta_records['errors'].append(error_dict)

                    
                # check that its compartment
                parent_node = morph.node_by_id(parent_id)
                parent_type = parent_node['type']

                if parent_type != cur_node_type:
                    # This is okay if it's parent is the soma
                    if parent_type != 1:
                        # Otherwise, this is unacceptable UNLESS it's axon stem from basasl
                        if (cur_node_type == 2) & (parent_type == 3):
                            axon_origins += 1
                            axon_origin_node_list.append(no)
                        else:
                            error_list.append(
                                "Node type {} has parent node of type {}".format(cur_node_type, parent_type))
                            node_list.append(no)
                            error_dict = {
                                'error':"Node type {} has parent node of type {}".format(cur_node_type, parent_type), 
                                  "Nodes":[no]
                                  }
                            meta_records['errors'].append(error_dict)


                    else:
                        # This node is a child of the soma
                        if cur_node_type == 2:
                            axon_origins += 1
                            axon_origin_node_list.append(no)

                        # Make sure it's not a furcation node
                        cur_node_children = morph.get_children(no)
                        if len(cur_node_children) > 1:
                            error_list.append("Node {} is an immediate child of the soma and branches".format(no['id']))
                            node_list.append(no)
                            error_dict = {
                                'error':"Node {} is an immediate child of the soma and branches".format(no['id']), 
                                  "Nodes":[no]
                                  }
                            meta_records['errors'].append(error_dict)

                            
                        # And make sure it's not too far away from the soma
                        if soma:

                            dist = euclidean(cur_node_coord, soma_coord)
                            if dist > distance_threshold:
                                error_list.append("Soma child node {} is {} distance from soma".format(no['id'], dist))
                                node_list.append(no)
                                error_dict = {
                                'error':"Soma child node {} is {} distance from soma".format(no['id'], dist), 
                                  "Nodes":[no]
                                  }
                                meta_records['errors'].append(error_dict)


            else:
                error_list.append(
                    "Node {}, Type {}, Parent ID {} Not In Morphology".format(no['id'], no['type'], parent_id))
                node_list.append(no)
                
                error_dict = {
                    'error':"Node {}, Type {}, Parent ID {} Not In Morphology".format(no['id'], no['type'], parent_id), 
                        "Nodes":[no]
                        }
                meta_records['errors'].append(error_dict)
                
        else:
            number_of_roots += 1
            root_node_list.append(no)

    if axon_origins > 1:
        error_list.append("Multiple Axon Origins ({} found)".format(axon_origins))
        [node_list.append(i) for i in axon_origin_node_list]
        error_dict = {
            'error':"Multiple Axon Origins ({} found)".format(axon_origins), 
            "Nodes":axon_origin_node_list
                }
        meta_records['errors'].append(error_dict)

        
    if number_of_roots > 1:
        error_list.append("Multiple Root Nodes ({} found)".format(number_of_roots))
        [node_list.append(i) for i in root_node_list]
        
        error_dict = {
            'error':"Multiple Root Nodes ({} found)".format(number_of_roots), 
            "Nodes":root_node_list
                }
        meta_records['errors'].append(error_dict)


    duplicate_coords = [k for k, v in number_of_nodes_at_each_coord.items() if v != 1]
    duplicate_nodes = [morph.node_by_id(node_coord_to_id[c]) for c in duplicate_coords]
    num_duplicate_coords = len(duplicate_coords)
    if num_duplicate_coords != 0:
        error_list.append("Nodes With Identical X,Y,Z Coordinates Found ({} found)".format(num_duplicate_coords))
        error_dict = {
            'error':"Nodes With Identical X,Y,Z Coordinates Found ({} found)".format(num_duplicate_coords), 
            "Nodes":duplicate_nodes
                }
        meta_records['errors'].append(error_dict)
        
        
        for dup_coord in duplicate_coords:
            dup_no_id = node_coord_to_id[dup_coord]
            this_no = morph.node_by_id(dup_no_id)
            node_list.append(this_no)
            
    # Loop Check
    has_loops, confidence = check_for_loops(morph)
    if has_loops:
        error_list.append("Loop Found In Morphology")
        error_dict = {
            'error':"Loop Found In Morphology", 
            "Nodes":[]
                }
        meta_records['errors'].append(error_dict)

    else:
        if confidence == "Ambiguous":
            error_list.append("Unable to check for loops due to missing root")
            error_dict = {
            'error':"Unable to check for loops due to missing root", 
            "Nodes":[]
                }
            meta_records['errors'].append(error_dict)

    # Soma ID Check:
    if (soma is not None) and (soma['id'] != 1):
        error_list.append("Soma Node ID Is Not 1")
        node_list.append(soma)
        error_dict = {
            'error':"Soma ID Is Not 1", 
            "Nodes":[soma]
                }
        meta_records['errors'].append(error_dict)

    
    meta_records['file_name'] = os.path.abspath(input_swc)
    
    error_list = list(set(error_list))

    res_dict = {"file_name": os.path.abspath(input_swc),
                "error_list": error_list}
    return meta_records


def check_for_loops(morphology):
    """
    Will use depth first traversal to identify any loops. If the morphology is not uniquely rooted, will visit
    each root node and check for loops in that tree.

    :param morphology:
    :return: loop status (True-has loops, False-none found), confidence (if ambiguous-unable to
    check for loops due to no clear starting point); Bool,str.
    """
    soma = morphology.get_soma()
    if not soma:
        start_nodes = [n for n in morphology.nodes() if n['parent'] == -1]
        if start_nodes == []:
            # Starting point ambiguous, may not be able to find loop
            print("WARNING: Unable to identify root node(s) for loop checking. Returning False")
            return False, "Ambiguous"
    else:
        start_nodes = [soma]

    for st_node in start_nodes:

        loops_below = dfs_loop_check(morphology, st_node)
        if loops_below:
            return True, "Confident"

    return False, "Confident"


def multiple_soma_nodes(morph):
    """
    Will check for multiple nodes of type 1
    :param morph: neuron_morphology morphology object
    :return: error_list
    """
    soma = morph.get_soma()
    type_1_nodes = [n for n in morph.nodes() if (n['type'] == 1) and (n['id'] != soma['id'])]

    if type_1_nodes != []:
        return ['Multiple Nodes ({}) of Type 1'.format(len(type_1_nodes))]
    else:
        return []


def duplicate_node_qc(morph):
    """
    Will check if multiple nodes occupy the same coordinate
    :param morph: neuron_morphology morphology object
    :return: error_list;list
    """
    record_dict = defaultdict(int)
    for no in morph.nodes():
        no_coord = (no['x'], no['y'], no['z'])
        record_dict[no_coord] += 1

    unique_counts = {k: v for k, v in record_dict.items() if v != 1}
    error_list = []
    for coord, ct in unique_counts.items():
        error_list.append("Coord: {}, {} instances".format(coord, ct))

    return error_list


def node_type_qc(morph, expected_types=[1, 2, 3, 4]):
    """
    Will make sure all nodes are within the constrained expected types

    :param morph: neuron_morphology object
    :param expected_types: list of expected node types in the morphology
    :return: error_list;list
    """
    error_list = []
    for no in morph.nodes():
        no_type = no['type']
        if no_type not in expected_types:
            error_list.append('Unexpected node type {}'.format(no_type))

    return set(error_list)


def soma_children_qc(morph, distance_threshold=50):
    """
    Will make sure the soma's children are within a certain threshold from the soma and do not branch.
    :param morph: neuron_morphology Morphology object
    :param distance_threshold: maximum distance a valid child may be from soma
    :return: error_list;list
    """
    error_list = []
    soma = morph.get_soma()
    soma_coord = [soma['x'], soma['y'], soma['z']]
    children = morph.get_children(soma)
    for ch_no in children:
        ch_coord = [ch_no['x'], ch_no['y'], ch_no['z']]
        dist = euclidean(soma_coord, ch_coord)

        if dist > distance_threshold:
            error_list.append("Soma child {} is {} distance from soma".format(ch_no['id'], dist))

        grandchildren = morph.get_children(ch_no)
        if len(grandchildren) > 1:
            error_list.append("Somas child is a furcation node")

    return error_list


def morphology_parent_node_qc(morph, types_to_check=[2, 3, 4]):
    """
    Will make sure parent-child relationships are valid for nodes in an swc file. This will only visit node types found
    in the keys of accepted_dict. Where valid means the parent id is actually in the morphology, and the parent node's
    type is acceptable. For example an axon node (2) can have a parent of soma (1), axon (2), or basal dendrite (3) in
    inhibitory neurons.

    :param morph: neuron_morphology morphology object
    :param types_to_check: list; compartment types to check
    :return: error_list;list
    """

    nodes_to_qc = [n for n in morph.nodes() if n['type'] in types_to_check]
    all_node_ids = [n['id'] for n in morph.nodes()]
    error_list = []
    axon_origin_nodes = []
    axon_origins = 0
    for no in nodes_to_qc:
        parent_id = no['parent']
        if parent_id in all_node_ids:
            cur_node_type = no['type']
            # check that its compartment
            parent_node = morph.node_by_id(parent_id)
            parent_type = parent_node['type']

            if parent_type != cur_node_type:
                # This is okay if it's parent is the soma
                if parent_type != 1:
                    # Otherwise, this is unacceptable UNLESS it's axon stem from basasl
                    if (cur_node_type == 2) & (parent_type == 3):
                        axon_origins += 1
                        axon_origin_nodes.append(no)
                        
                    else:
                        error_list.append("Node type {} has parent node of type {}".format(cur_node_type, parent_type))

                else:
                    if cur_node_type == 2:
                        axon_origins += 1
            # acceptable_types = accepted_dict[no['type']]
            # if parent_type not in acceptable_types:
            #     error_list.append("Node type {} has parent node of type {}".format(no['type'], parent_type))

        else:
            error_list.append(
                "Node {}, Type {}, Parent ID {} Not In Morphology".format(no['id'], no['type'], parent_id))

    if axon_origins > 1:
        error_list.append("Multiple Axon Origins ({} found".format(axon_origins))

    error_list = list(set(error_list))
    return error_list
