import json
import numpy as np
from importlib.resources import files

def load_layer_template(species):
    """ Load the average cortical layer depth JSON file for a species.

    Keys are strings representing cortical layers (e.g.'2/3','4'...)
    Values represent the cortical depth for the top (pia side) of a given layer

    TODO add NHP average layers when generated. 

    :param species: species of average layers to return. Options: 'mouse' or 'human' 
    :return: layers, dictionary of distances to the pia (in microns) from the upper side of each layer.
    :return: labels, dictionary of layer names (e.g. 'L1', 'L2/3') and midpoints 
    :return: colors, dictionary of layer line colors 
    """
    
    #average layer depths 
    depth_file = files('morph_utils') / 'data/{}_average_layer_depths.json'.format(species)
    with open(depth_file, "r") as fn: 
        layers = json.load(fn)

    #layer line colors 
    layer_keys = list(layers.keys())
    layer_colors = ['darkgrey' if l in ['1', 'wm'] else 'lightgrey' for l in layer_keys]
    colors = dict(zip(layer_keys, layer_colors))

    #layer labels and midpoints 
    layer_depths = np.array(list(layers.values()))
    layer_midpoints = list((layer_depths[1:] + layer_depths[:-1]) / 2)
    layer_names = ['L' + l for l in layer_keys[:-1]]
    labels = dict(zip(layer_names, layer_midpoints))

    layer_dict = {'layers': layers, 
                  'labels': labels, 
                  'colors': colors}

    return layer_dict 