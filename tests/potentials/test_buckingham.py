import unittest
import pyflamestk.potentials as potentials

class TestBuckinghamInitialize(unittest.TestCase):
    def setUp(self):
        self.name            = 'buckingham test'
        self.atoms           = ['Mg','O']
        self.parameter_names = ['chrg_Mg', 'chrg_O', 
                                'MgMg_A', 'MgMg_rho', 'MgMg_C', 
                                'MgO_A', 'MgO_rho', 'MgO_C', 
                                'OO_A', 'OO_rho', 'OO_C']
        self.parameter_constraints = [['MgMg_rho', '>=', 0],
                                      ['MgO_rho', '>=', 0],
                                      ['OO_rho', '>=', 0]]         
        self.potential = potentials.Buckingham(self.name,
                                               self.atoms)

    def test_parameter_names(self):
        self.assertEqual(self.potential.parameter_names,
                         self.parameter_names)
                         
    def test_parameter_dictionary(self):
        dict_params = {p for p in self.potential._params.keys()}
        obj_params  = {p for p in self.potential.parameter_names}
        self.assertEqual(dict_params,obj_params)
    
    def test_parameter_constraints(self):
        self.assertEqual(len(self.parameter_constraints),
                         len(self.potential.parameter_constraints))
        for c in self.potential.parameter_constraints:
            self.assertIn(c,self.parameter_constraints)
        
                          
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBuckinghamInitialize)
    unittest.TextTestRunner(verbosity=2).run(suite)