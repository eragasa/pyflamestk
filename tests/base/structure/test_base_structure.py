import copy
import unittest
import numpy as np
import numpy.testing as npt
import pyflamestk.base as base

class TestStructureInitialize(unittest.TestCase):
    def setUp(self):
        self.a0 = 1.0
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
        self.uc_1.lattice_parameter = self.a0
        self.uc_1.h_matrix = copy.deepcopy(self.h_matrix)
        for a in self.atoms:
            self.uc_1.add_atom(a[0],copy.deepcopy(a[1]))

    def test_lattice_parameter(self):
        self.assertEqual(self.uc_1.lattice_parameter, self.a0)
        
    def test_h_matrix_is_equal(self):
        npt.assert_allclose(self.h_matrix,
                            self.uc_1.h_matrix)

    def test_n_atoms(self):
        n_atoms = len(self.atoms)
        self.assertEqual(n_atoms,self.uc_1.n_atoms)

    def test_n_iterstitials(self):
        n_iterstitials = 0
        self.assertEqual(n_iterstitials,len(self.uc_1.interstitials))

    def test_n_vacancies(self):
        n_vacancies = 0
        self.assertEqual(n_vacancies,len(self.uc_1.vacancies))

    def test_h1(self):
        H = np.array(self.h_matrix)
        h1 = H[0,:]
        npt.assert_allclose(h1,self.uc_1.h1)

    def test_h1(self):
        H = np.array(self.h_matrix)
        h2 = H[1,:]
        npt.assert_allclose(h2,self.uc_1.h2)

    def test_h1(self):
        H = np.array(self.h_matrix)
        h3 = H[2,:]
        npt.assert_allclose(h3,self.uc_1.h3)

    def test_a0(self):
        self.assertEqual(self.a0,self.uc_1.lattice_parameter)

    def test_a1(self):
        H = np.array(self.h_matrix)
        a0 = H[0,0]
        npt.assert_allclose(a0,self.uc_1.a1)

    def test_a2(self):
        H = np.array(self.h_matrix)
        a0 = H[0,0]
        self.assertEqual(a0,self.uc_1.a2)

    def test_a3(self):
        H = np.array(self.h_matrix)
        a0 = H[0,0]
        self.assertEqual(a0,self.uc_1.a3)

    def TearDown(): pass    
            

class Test_Base_Structure__AddAtom(unittest.TestCase):
    def setUp(self):
        self.a0 = 1.0
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
        self.uc_1.lattice_parameter = self.a0
        self.uc_1.h_matrix = copy.deepcopy(self.h_matrix)
        for a in self.atoms:
            self.uc_1.add_atom(a[0],copy.deepcopy(a[1]))

    def test_add_atom_to_available_position(self):
        self.uc_1.add_atom('Mg',[1/4,1/4,1/4])
        n_atoms = len(self.atoms) + 1
        self.assertEqual(n_atoms,self.uc_1.n_atoms)

    def test_add_atom_available_position_check_n_vacancies_are_zero(self):
        self.uc_1.add_atom('Mg',[1/4,1/4,1/4])
        self.assertEqual(0,len(self.uc_1.vacancies))

    def test_add_atom_available_position_check_n_interstitials_are_zero(self):
        self.uc_1.add_atom('Mg',[1/4,1/4,1/4])
        self.assertEqual(0,len(self.uc_1.interstitials))

    def test_add_atom_to_unavailable_position_raises_error(self):
        self.assertRaises(ValueError,
                          lambda: self.uc_1.add_atom('Mg',[0,0,0]))

    def test_add_atom_to_unavailable_position_check_n_atoms(self):
        try:
            self.uc_1.add_atom('Mg',[0,0,0])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(n_atoms,self.uc_1.n_atoms)
    def test_add_atom_available_position_check_n_vacancies_are_zero(self):
        try:
            self.uc_1.add_atom('Mg',[0,0,0])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(0,len(self.uc_1.vacancies))

    def test_add_atom_available_position_check_n_interstitials_are_zero(self):
        try:
            self.uc_1.add_atom('Mg',[0,0,0])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(0,len(self.uc_1.interstitials))

class TestBaseStructureRemoveAtom(unittest.TestCase):
    def setUp(self):
        self.a0 = 1.0
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
        self.uc_1.lattice_parameter = self.a0
        self.uc_1.h_matrix = copy.deepcopy(self.h_matrix)
        for a in self.atoms:
            self.uc_1.add_atom(a[0],copy.deepcopy(a[1]))

    def test_remove_atom_from_unoccupied_position_raises_error(self):
        self.assertRaises(ValueError,
                          lambda: self.uc_1.remove_atom('Mg',[1/4,1/4,1/4]))

    def test_remove_atom_to_unoccupied_position_check_n_atoms(self):
        try:
            self.uc_1.remove_atom('Mg',[1/4,1/4,1/4])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(n_atoms,self.uc_1.n_atoms)

    def test_remove_atom_unoccupied_position_check_n_vacancies_are_zero(self):
        try:
            self.uc_1.remove_atom('Mg',[1/4,1/4,1/4])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(0,len(self.uc_1.vacancies))

    def test_remove_atom_unoccupied_position_check_n_interstitials_are_zero(self):
        try:
            self.uc_1.remove_atom('Mg',[1/4,1/4,1/4])
        except ValueError:
            pass
        n_atoms = len(self.atoms)
        self.assertEqual(0,len(self.uc_1.interstitials))

    def test_remove_atom_to_occupied_position_check_n_atoms(self):
        self.uc_1.remove_atom('Mg',[1/2,1/2,0])
        n_atoms = len(self.atoms) - 1
        self.assertEqual(n_atoms,self.uc_1.n_atoms)

    def test_remove_atom_occupied_position_check_n_vacancies_are_zero(self):
        self.uc_1.remove_atom('Mg',[1/2,1/2,0])
        self.assertEqual(0,len(self.uc_1.vacancies))

    def test_remove_atom_available_position_check_n_interstitials_are_zero(self):
        self.uc_1.remove_atom('Mg',[1/2,1/2,0])
        self.assertEqual(0,len(self.uc_1.interstitials))

if __name__ == '__main__':
    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructureInitialize)
    unittest.TextTestRunner(verbosity=2).run(suite)
