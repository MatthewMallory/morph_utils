from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
from sklearn.neighbors import KDTree
from neuron_morphology.swc_io import morphology_from_swc

def cellwidth(morph, compartments=[2, 3, 4]):
    right_extent = rightextent(morph, compartments)
    left_extent = leftextent(morph, compartments)
    return right_extent + left_extent


def downextent(morph, compartments=[2, 3, 4]):
    return min([n['y'] - morph.get_soma()['y'] for n in morph.nodes() if n['type'] in compartments] + [0])


def upextent(morph, compartments=[2, 3, 4]):
    return max([n['y'] - morph.get_soma()['y'] for n in morph.nodes() if n['type'] in compartments] + [0])


def leftextent(morph, compartments=[2, 3, 4]):
    return abs(min([n['x'] - morph.get_soma()['x'] for n in morph.nodes() if n['type'] in compartments] + [0]))


def rightextent(morph, compartments=[2, 3, 4]):
    return max([n['x'] - morph.get_soma()['x'] for n in morph.nodes() if n['type'] in compartments] + [0])


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


def get_node_spacing(morph, node_types_to_check=[1, 3, 4]):
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


def node_distance_between_morphs(swc_file_1, swc_file_2, compartment_types):
    """
    will calculate the mean distance between two swc file nodes of a certain type. For each  node in file 1, the
    nearest node in file 2 is found (that has type contained in compartment_types). This is the forward distance. The
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
