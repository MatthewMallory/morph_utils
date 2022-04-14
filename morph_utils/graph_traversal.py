from collections import deque

def bfs_tree(st_node, morph):
    """
    breadth first traversal of tree, returns nodes in segment and how many

    :param st_node: node to begin BFS traversal from.
    :param morph: neuron_morphology Morphology object
    :return: list of nodes in segment (including start node), int number of nodes in segment
    """
    queue = deque([st_node])
    nodes_in_segment = []
    seg_len = 0
    while len(queue) > 0:
        seg_len += 1
        current_node = queue.popleft()
        nodes_in_segment.append(current_node)
        for ch_no in morph.get_children(current_node):
            queue.append(ch_no)

    return nodes_in_segment, len(nodes_in_segment)

def dfs_labeling(st_node, new_starting_id, modifying_dict, morph):
    """
    depth first traversal for relabeling a segment of a morphology.
    :param st_node: a node to start labelling from
    :param new_starting_id: the new node id to assign said start node
    :param modifying_dict: retains information on node id updates
    :param morph: neuron_morphology Moprhology object
    :return:
    """
    ct = 0
    queue = deque([st_node])
    while len(queue) > 0:
        ct += 1
        current_node = queue.popleft()
        modifying_dict[current_node['id']] = new_starting_id
        new_starting_id += 1
        for ch_no in morph.get_children(current_node):
            queue.appendleft(ch_no)
    return ct

def get_path_to_root(start_node, morphology):
    """
    get the nodes along the path from a given start node to a root node (where root node has parent = -1)
    :param start_node:
    :param morphology:
    :return: list of nodes including start node
    """
    seg_up = []
    current_node = start_node
    seg_up.append(start_node)
    current_parent_id = current_node['parent']
    iteration_count = 0
    max_iterations = len(morphology.nodes())
    while current_parent_id != -1:
        current_node = morphology.node_by_id(current_parent_id)
        seg_up.append(current_node)
        current_parent_id = current_node['parent']
        iteration_count+=1

        if iteration_count > max_iterations:
            print("Iterations Exceeded number of nodes in morphology. Check for loops with "
                  "morph_utils.traversal.check_for_loops(morphology)")
            return None
    return seg_up




def check_for_loops(morphology):
    #check for roots
    return True
