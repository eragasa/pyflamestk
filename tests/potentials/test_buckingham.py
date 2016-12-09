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
        # test that you have the correct number of constraints
        self.assertEqual(len(self.parameter_constraints),
                         len(self.potential.parameter_constraints))
                         
        # test that there are no double constraints
        self.assertCountEqual(self.parameter_constraints,
                              self.potential.parameter_constraints,
                              1)
                                  
        # test the each constraint is supposed to be there
        for c in self.potential.parameter_constraints:
            self.assertIn(c,self.parameter_constraints)

class TestBuckinghamAddConstraint(unittest.TestCase):
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
                                               
    def test_add_same_constraint(self):
        pass

    def test_add_different_constraint(self):                                               
        new_c = ['chrg_Mg', '=', '-chrg_O']
        self.potential.add_parameter_constraint(new_c[0],new_c[1],new_c[2])
        # test that there are no double constraints
        self.assertCountEqual(self.parameter_constraints,
                              self.potential.parameter_constraints,
                              1)

        
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBuckinghamInitialize)
    suite.addTest(TestBuckinghamAddConstraint)
    unittest.TextTestRunner(verbosity=2).run(suite)