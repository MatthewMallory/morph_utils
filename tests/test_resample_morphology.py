import os
import unittest
from neuron_morphology.swc_io import Morphology
from morph_utils.modifications import resample_morphology
from morph_utils.measurements import get_node_spacing

class TestResample(unittest.TestCase):

    def setUp(self):
        dummy_nodes =[
                {
                "x":0,"y":0,"z":0,
                'id':1,
                "type":1,"radius":1,
                "parent":-1,
                },
                {
                "x":0,"y":-10,"z":0,
                'id':2,
                "type":2,"radius":1,
                "parent":1,
                },
                {
                "x":0,"y":-15,"z":0,
                'id':3,
                "type":2,"radius":1,
                "parent":2,
                },
                {
                "x":0,"y":-25,"z":0,
                'id':4,
                "type":2,"radius":1,
                "parent":3,
                },
                
                
                {
                "x":0,"y":-30,"z":0,
                'id':5,
                "type":2,"radius":1,
                "parent":4,
                },
                {
                "x":15,"y":-35,"z":0,
                'id':6,
                "type":2,"radius":1,
                "parent":5,
                },
                {
                "x":20,"y":-45,"z":0,
                'id':7,
                "type":2,"radius":1,
                "parent":6,
                },
                
                {
                "x":-15,"y":-35,"z":0,
                'id':8,
                "type":2,"radius":1,
                "parent":4,
                },
                {
                "x":-55,"y":-35,"z":0,
                'id':9,
                "type":2,"radius":1,
                "parent":8,
                },
                
                
                {
                "x":0,"y":20,"z":0,
                'id':10,
                "type":2,"radius":1,
                "parent":1,
                },
                {
                "x":0,"y":45,"z":0,
                'id':11,
                "type":2,"radius":1,
                "parent":10,
                },
                {
                "x":25,"y":55,"z":0,
                'id':12,
                "type":2,"radius":1,
                "parent":11,
                },
                {
                "x":-25,"y":55,"z":0,
                'id':13,
                "type":2,"radius":1,
                "parent":11,
                },
        ]

        dumy_morph = Morphology(dummy_nodes,
                                parent_id_cb=lambda x:x['parent'],
                                node_id_cb=lambda x:x['id'])

        self.morphology = dumy_morph
        
    def test_upsampling(self):
        
        # test upsampling
        morph = self.morphology.clone()
        upsampled_morph = resample_morphology(morph, 0.25)
        up_spacing = get_node_spacing(upsampled_morph)[0]
        
        # test downsampling
        morph = self.morphology.clone()
        downsampled_morph = resample_morphology(morph, 30)    
        down_spacing = get_node_spacing(downsampled_morph)[0]
        
        self.assertAlmostEqual(up_spacing,0.28938159542643754)
        self.assertAlmostEqual(down_spacing,20.90397314689175)


if __name__ == '__main__':
    unittest.main()
            