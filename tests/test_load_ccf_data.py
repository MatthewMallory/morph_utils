import unittest
from morph_utils.ccf import open_ccf_annotation, load_structure_graph
import SimpleITK as sitk
import numpy as np

class TestImportFunctions(unittest.TestCase):

    def test_imports(self):

        sg_df = load_structure_graph()
        self.assertTrue(not sg_df.empty)
        # annotation_array = open_ccf_annotation(as_array=True)
        # self.assertTrue(isinstance(annotation_array, np.ndarray))

if __name__ == '__main__':
    unittest.main()
