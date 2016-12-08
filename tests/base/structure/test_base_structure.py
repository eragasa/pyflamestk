import copy
import unittest
import numpy as np
import pyflamestk.base as base

class TestStructureInitialize(unittest.TestCase):
    def setUp(self):
        self.a = 1.0
        self.h_matrix = [[4.256483, 0.000000, 0.000000],
                         [0.000000, 4.256483, 0.000000],
                         [0.000000, 0.000000, 4.256483]]
        self.atoms = [['Mg',[0.500000, 0.500000, 0.000000]],
                      ['Mg',[0.500000, 0.000000, 0.500000]], 
                      ['Mg',[0.000000, 0.000000, 0.000000]],
                      ['Mg',[0.000000, 0.500000, 0.500000]],
                      ['O', [0.000000, 0.500000, 0.000000]],
                      ['O', [0.500000, 0.000000, 0.000000]],
                      ['O', [0.500000, 0.500000, 0.500000]],
                      ['O', [0.000000, 0.000000, 0.500000]]]

        self.uc_1 = base.Structure()
        self.uc_1.lattice_parameter = self.a
        self.uc_1.h_matrix = copy.deepcopy(self.h_matrix)
        for a in self.atoms:
            self.uc_1.add_atom(a[0],copy.deepcopy(a[1]))

    def test_lattice_parameter(self):
        self.assertEqual(self.uc_1.lattice_parameter, self.a)
        
    def test_h_matrix(self):
        self.assertTrue(np.testing.assert_array_equal(self.uc_1.h_matrix,
                                                      self.h_matrix))
        
    def TearDown(): pass    
            

#class TestStructureCopy(unittest.TestCase):
#    
#    def setUp():
#        uc_2 = base.Structure()
#        uc_2.use_structure(uc_1)
#
#    
#    def tearDown():
#        pass
    
if __name__ == '__main__':
    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructureInitialize)
    unittest.TextTestRunner(verbosity=2).run(suite)