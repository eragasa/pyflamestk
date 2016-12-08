import unittest

class TestOutcar(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(),'FOO')
        
    
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestOutcar)
    unittest.TextTestRunner(verbosity=2).run(suite)