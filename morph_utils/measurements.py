from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
from sklearn.neighbors import KDTree
from neuron_morphology.swc_io import morphology_from_swc


def cellwidth(morph, compartments=[2, 3, 4]):
    right_extent = rightextent(morph, compartments)
    left_extent = leftextent(morph, compartments)
    return right_extent + left_extent


def cellheight(morph, compartments=[2, 3, 4]):
    up = upextent(morph, compartments)
    down = downextent(morph, compartments)
    return up + down


def downextent(morph, compartments=[2, 3, 4]):
    return min([n['y'] - morph.get_soma()['y'] for n in morph.nodes() if n['type'] in compartments] + [0])


def upextent(morph, compartments=[2, 3, 4]):
    return max([n['y'] - morph.get_soma()['y'] for n in morph.nodes() if n['type'] in compartments] + [0])


def leftextent(morph, compartments=[2, 3, 4]):
    return abs(min([n['x'] - morph.get_soma()['x'] for n in morph.nodes() if n['type'] in compartments] + [0]))


def rightextent(morph, compartments=[2, 3, 4]):
    return max([n['x'] - morph.get_soma()['x'] for n in morph.nodes() if n['type'] in compartments] + [0])


def tree_length(morphology, st_node=None):
    "Measure the total length of the descendant tree from a given start node"

    if not st_node:
        st_node = morphology.get_soma()

        if not st_node:
            st_node = [n for n in morphology.nodes() if n['parent'] == -1][0]

    visited = []
    dist = 0.0
    queue = deque([st_node])
    while len(queue) > 0:
        current_node = queue.popleft()
        visited.append(current_node)
        for ch_no in morphology.get_children(current_node)[::-1]:
            queue.appendleft(ch_no)
            dist += dist_bwn_nodes(current_node, ch_no)
    return dist


def dist_bwn_nodes(n1, n2):
    """
    Euclidean distance between nodes
    :param n1: node
    :param n2: node
    :return: float, average euclidean distance
    """
    try:
        return euclidean((n1['x'], n1['y'], n1['z']),
                         (n2['x'], n2['y'], n2['z']))
    except:
        return np.nan


def get_node_spacing(morph, node_types_to_check=[1, 2, 3, 4]):
    """
    Will calculate the average spacing between nodes in a morphology.

    :param morph: neuron_morphology morphology
    :param node_types_to_check: list of ints, which compartment(s) to consider
    :return: (float, list) average spacing and list of the raw values
    """
    st_node = morph.get_soma()

    queue = deque([st_node])
    nodes_in_segment = []
    distances = []
    while len(queue) > 0:
        current_node = queue.popleft()
        nodes_in_segment.append(current_node)
        for ch_no in morph.get_children(current_node):
            if ch_no['type'] in node_types_to_check:
                dist_btwn = dist_bwn_nodes(current_node, ch_no)
                distances.append(dist_btwn)
                queue.append(ch_no)

    return np.nanmean(distances), distances


def simple_node_distance_between_morphs(swc_file_1, swc_file_2, compartment_types):
    """
    will calculate the mean distance between two swc file nodes of a certain type. For each  node in file 1, the
    nearest node in file 2 is found (that has type contained in compartment_types). The matched node does not need
    to be of the same compartment, only needs to be found in compartment_types. This is the forward distance. The
    reverse distance will be the opposite direction. Return is a dictionary for easy multiprocessing compatibility

    Note that basal dendrites can be matched to apical dendrites and vice versa if compartment_types = [3,4] if that
    is the closest coordinate found.

    :param swc_file_1: str, path to swc file 1
    :param swc_file_2: str, path to swc file 2
    :param compartment_types: list of ints, which compartments you want used in constructing distance measurements
    :return: dict, dictionary with forward (file 1 -> file 2) and reverse distances (file 2 -> file 1)
    """
    morph_1 = morphology_from_swc(swc_file_1)
    morph_2 = morphology_from_swc(swc_file_2)

    morph_1_nodes = np.array([[n['x'], n['y'], n['z']] for n in morph_1.nodes() if n['type'] in compartment_types])
    morph_2_nodes = np.array([[n['x'], n['y'], n['z']] for n in morph_2.nodes() if n['type'] in compartment_types])

    morph_1_kd_tree = KDTree(morph_1_nodes)
    morph_2_kd_tree = KDTree(morph_2_nodes)

    morph_2_to_1_dists, _ = morph_1_kd_tree.query(morph_2_nodes, k=1)
    morph_1_to_2_dists, _ = morph_2_kd_tree.query(morph_1_nodes, k=1)

    forward_mean_dist = np.mean(morph_1_to_2_dists)
    reverse_mean_dist = np.mean(morph_2_to_1_dists)

    results_dict = {"file_1": swc_file_1,
                    "file_2": swc_file_2,
                    "forward_distance": forward_mean_dist,
                    "reverse_distance": reverse_mean_dist}
    return results_dict


def node_distance_between_morphs(swc_file_1, swc_file_2, compartment_types, compartment_match=True):
    """
    will calculate the mean distance between two swc file nodes of a certain type. For each  node in file 1, the
    nearest node in file 2 is found (that has type contained in compartment_types). When compartment_match is True, the
    matched node from file 2 must be of the same compartment as the query node in file 1, otherwise the score is
    penalized This is the forward distance. The reverse distance will be the opposite direction. Return is a dictionary
    for easy multiprocessing compatibility

    Note that basal dendrites can be matched to apical dendrites and vice versa if compartment_types = [3,4] if that
    is the closest coordinate found AND compartment_match = False. When compartment_match is True, a penalty is added
    to the distance score if compartments are mismatched.

    :param swc_file_1: str, path to swc file 1
    :param swc_file_2: str, path to swc file 2
    :param compartment_types: list of ints, which compartments you want used in constructing distance measurements
    :param compartment_match: bool, if True, will only match dendrite to dendrite, axon to axon, and soma to soma nodes
    :return: dict, dictionary with forward (file 1 -> file 2) and reverse distances (file 2 -> file 1)
    """
    morph_1 = morphology_from_swc(swc_file_1)
    morph_2 = morphology_from_swc(swc_file_2)

    all_morph_1_nodes = [n for n in morph_1.nodes() if n['type'] in compartment_types]
    all_morph_2_nodes = [n for n in morph_2.nodes() if n['type'] in compartment_types]

    if len(all_morph_1_nodes) == 0 or len(all_morph_2_nodes) == 0:
        return {"file_1": swc_file_1,
                "file_2": swc_file_2,
                "forward_distance": np.inf,
                "reverse_distance": np.inf}

    all_morph_1_nodes_arr = np.array([[n['x'], n['y'], n['z']] for n in all_morph_1_nodes])
    all_morph_2_nodes_arr = np.array([[n['x'], n['y'], n['z']] for n in all_morph_2_nodes])

    if compartment_match:

        morph_2_to_1_all_dists, morph_1_to_2_all_dists = [], []
        for comp_type in compartment_types:
            morph_1_nodes = np.array([[n['x'], n['y'], n['z']] for n in morph_1.nodes() if n['type'] == comp_type])
            morph_2_nodes = np.array([[n['x'], n['y'], n['z']] for n in morph_2.nodes() if n['type'] == comp_type])

            if len(morph_1_nodes) > 0 and len(morph_2_nodes) > 0:
                morph_1_kd_tree = KDTree(morph_1_nodes)
                morph_2_kd_tree = KDTree(morph_2_nodes)

                morph_2_to_1_dists, _ = morph_1_kd_tree.query(morph_2_nodes, k=1)
                morph_1_to_2_dists, _ = morph_2_kd_tree.query(morph_1_nodes, k=1)

                [morph_2_to_1_all_dists.append(d) for d in morph_2_to_1_dists]
                [morph_1_to_2_all_dists.append(d) for d in morph_1_to_2_dists]

            elif len(morph_1_nodes) == 0 and len(morph_2_nodes) == 0:
                continue
            else:

                # find nearest non compartment match and penalize
                if len(morph_1_nodes) == 0:

                    morph_1_kd_tree = KDTree(all_morph_1_nodes_arr)
                    morph_2_kd_tree = KDTree(morph_2_nodes)

                    morph_2_to_1_dists, _ = morph_1_kd_tree.query(morph_2_nodes, k=1)
                    morph_1_to_2_dists, _ = morph_2_kd_tree.query(all_morph_1_nodes_arr, k=1)

                else:
                    morph_1_kd_tree = KDTree(morph_1_nodes)
                    morph_2_kd_tree = KDTree(all_morph_2_nodes_arr)

                    morph_2_to_1_dists, _ = morph_1_kd_tree.query(all_morph_2_nodes_arr, k=1)
                    morph_1_to_2_dists, _ = morph_2_kd_tree.query(morph_1_nodes, k=1)

                # penalize distance metric for incompatible compartment matching
                [morph_2_to_1_all_dists.append(d * 3.0) for d in morph_2_to_1_dists]
                [morph_1_to_2_all_dists.append(d * 3.0) for d in morph_1_to_2_dists]

    else:

        morph_1_kd_tree = KDTree(all_morph_1_nodes_arr)
        morph_2_kd_tree = KDTree(all_morph_2_nodes_arr)

        morph_2_to_1_all_dists, _ = morph_1_kd_tree.query(all_morph_2_nodes_arr, k=1)
        morph_1_to_2_all_dists, _ = morph_2_kd_tree.query(all_morph_1_nodes_arr, k=1)

    forward_mean_dist = np.mean(morph_1_to_2_all_dists)
    reverse_mean_dist = np.mean(morph_2_to_1_all_dists)

    results_dict = {"file_1": swc_file_1,
                    "file_2": swc_file_2,
                    "forward_distance": forward_mean_dist,
                    "reverse_distance": reverse_mean_dist}
    return results_dict
