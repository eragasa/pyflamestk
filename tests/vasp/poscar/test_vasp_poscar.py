import unittest
import pyflamestk.base as base
import pyllamestk.vasp as vasp

class TestPoscarCopyMethod(unittest.TestCase):
    
    def setUp():
        uc = vasp.Poscar()
        
    def test_upper(self):
        
        self.assertEqual('foo'.upper(),'FOO')
    
    
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestOutcar)
    unittest.TextTestRunner(verbosity=2).run(suite)