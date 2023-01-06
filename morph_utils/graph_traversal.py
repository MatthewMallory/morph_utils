from collections import deque
import numpy as np
from scipy.spatial.distance import euclidean

def bfs_tree(st_node, morph):
    """
    breadth first traversal of tree, returns nodes in segment and how many.

    :param st_node: node to begin BFS traversal from.
    :param morph: neuron_morphology Morphology object
    :return: list of nodes in segment (including start node), int number of nodes in segment
    """
    max_iterations = len(morph.nodes())
    queue = deque([st_node])
    nodes_in_segment = []
    seg_len = 0
    while len(queue) > 0:
        seg_len += 1
        if seg_len > max_iterations:
            return [], 0
        current_node = queue.popleft()
        nodes_in_segment.append(current_node)
        for ch_no in morph.get_children(current_node):
            queue.append(ch_no)

    return nodes_in_segment, len(nodes_in_segment)


def dfs_tree(morphology, st_node):
    """
    Graph traversal using depth first search
    :param morphology: neuron_morphology.Morphology object
    :param st_node: neuron_morphology.Morphology.node, dictionary
    :return: list of nodes in segment (in dfs order including start node), int number of nodes in segment
    """
    visited = []
    queue = deque([st_node])
    while len(queue) > 0:
        current_node = queue.popleft()
        visited.append(current_node)
        for ch_no in morphology.get_children(current_node):
            queue.appendleft(ch_no)

    return visited, len(visited)


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




def dfs_loop_check(morphology, st_node):
    """
    Given a starting point on the graph, are there any loops in the subsequent tree.

    :param st_node: node to begin BFS traversal from; dict
    :param morphology: neuron_morphology Morphology object
    :return: True if there is a loop, False if not; Bool
    """
    visited = set()
    queue = deque([st_node])
    while len(queue) > 0:
        current_node = queue.popleft()
        current_node_id = current_node['id']
        if current_node_id in visited:
            return True
        else:
            visited.add(current_node_id)
            for ch_no in morphology.get_children(current_node):
                queue.appendleft(ch_no)

    return False


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


def get_path_and_path_dist_between_two_nodes(lower_node, upper_node, morphology):
    """
    Will return the path distance between two nodes assuming they are on the same branch. If they are not on the same
    branch, will return np.inf. In the future may refactor this so that we first check both nodes are in the tree,
    then if they are of different branches, return the distance from lower_node to common ancestor and upper_node to
    common ancestor. That may not fall within the scope of this function though.

    :param lower_node: neuron_morphology node that is a descendent of upper_node; dict
    :param upper_node: neuron_morphology node that is an ancestor to lower_node; dict
    :param morphology: neuron_morphology morphology
    :return: the path distance between the input nodes; float
    """
    upper_node_id = upper_node['id']
    lower_parent_id = lower_node['parent']

    max_iterations = len(morphology.nodes())
    path_distance = 0.0
    path = []
    counter = 0
    while lower_node['id'] != upper_node_id:
        counter += 1
        if counter > max_iterations:
            return [], np.inf

        parent_node = morphology.node_by_id(lower_parent_id)
        path.append(lower_node)
        path_distance += euclidean([parent_node['x'], parent_node['y'], parent_node['z']],
                                   [lower_node['x'], lower_node['y'], lower_node['z']])

        lower_node = parent_node
        lower_parent_id = lower_node['parent']

    return path, path_distance
