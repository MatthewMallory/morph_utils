import os
import unittest
from neuron_morphology.swc_io import Morphology, morphology_to_swc
from morph_utils.ccf import projection_matrix_for_swc, open_ccf_annotation

class TestProjectionMatrix(unittest.TestCase):

    def setUp(self):

        soma_x,soma_y,soma_z = 7992.894, 1786.669, 9253.204,
        morphology = Morphology([
                        {
                            "id": 0,
                            "parent": -1,
                            "type": 1,
                            "x": soma_x,
                            "y": soma_y,
                            "z": soma_z,
                            "radius": 10
                        },

                        {
                            "id": 1,
                            "parent": 0,
                            "type": 2,
                            "x": soma_x,
                            "y": soma_y,
                            "z": soma_z+10,
                            "radius": 3
                        },
                        {
                            "id": 2,
                            "parent": 1,
                            "type": 2,
                            "x": soma_x,
                            "y": soma_y+10,
                            "z": soma_z+10,
                            "radius": 3
                        },
                        {
                            "id": 3,
                            "parent": 1,
                            "type": 2,
                            "x": soma_x,
                            "y": soma_y-10,
                            "z": soma_z+10,
                            "radius": 3
                        },
                    ],
                    node_id_cb=lambda node: node["id"],
                    parent_id_cb=lambda node: node["parent"],
                )


        temp_ofile = "TempMorphFile.swc"
        morphology_to_swc(morphology, temp_ofile)

        self.swc_file = temp_ofile
        self.morphology = morphology
        self.annotation  = open_ccf_annotation(with_nrrd=True)
    
    def test_proj_mat_from_swc(self):
        swc_path = self.swc_file
        annotation = open_ccf_annotation(with_nrrd=True)
        projection_results = projection_matrix_for_swc(swc_path,
                                                        branch_count=False,
                                                        annotation=self.annotation , 
                                                        annotation_path = None, 
                                                        volume_shape=(1320, 800, 1140),
                                                        resolution=10)
        projection_dicts = list(projection_results)[1:]
        for d in projection_dicts:
            self.assertEqual(d['ipsi_VISal5'],30.0)

    def test_proj_mat_from_swc_branch_ct(self):
        swc_path = self.swc_file
        annotation = open_ccf_annotation(with_nrrd=True)
        projection_results = projection_matrix_for_swc(swc_path,
                                                        branch_count=True,
                                                        annotation=self.annotation , 
                                                        annotation_path = None, 
                                                        volume_shape=(1320, 800, 1140),
                                                        resolution=10)
        projection_dicts = list(projection_results)[1:]
        for d in projection_dicts:
            self.assertEqual(d['ipsi_VISal5'],30.0)

    
    def tearDown(self): 
        os.remove(self.swc_file)         
    
# if __name__ == '__main__':
#     unittest.main()
            