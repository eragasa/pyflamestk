import copy
import unittest
import numpy as np
import numpy.testing as npt
import collections # required to test dictionaries
import pyflamestk.base as base
import pyflamestk.pyposmat as pyposmat

class TestPyposmatQoiConfigFile(unittest.TestCase):
    def setUp(self):
        fname_qoi_config = 'pyposmat.qoi'
        self.qoi_config = pyposmat.QoiConfigFile(fname_qoi_config)
        self.qoi_targets = {'MgO_NaCl.a0':4.246,
                            'MgO_NaCl.c11':277.00,
                            'MgO_NaCl.c12':91.67,
                            'MgO_NaCl.c44':144.01,
                            'MgO_NaCl.B':153.45,
                            'MgO_NaCl.G':92.66,
                            'MgO_NaCl.fr_a':10.978,
                            'MgO_NaCl.fr_c':8.986,
                            'MgO_NaCl.sch':5.067,
                            'MgO_NaCl.001s':0.05595}
        self.qoi_names = ['MgO_NaCl.a0',
                          'MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
                          'MgO_NaCl.B', 'MgO_NaCl.G',
                          'MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch',
                          'MgO_NaCl.001s']
    def assertDictEqual(self, d1, d2, msg=None):
        # assertEqual uses for dicts
        for k,v1 in d1.items():
            self.assertIn(k, d2, msg)
            v2 = d2[k]
            if(isinstance(v1, collections.Iterable) and
                not isinstance(v1, basestring)):
                    self.assertItemsEqual(v1, v2, msg)
            else:
                self.assertEqual(v1, v2, msg)
                return True
    def test_attribute_qoi_target_is_equal(self):
        self.assertEqual(self.qoi_config.qoi_targets,
                         self.qoi_targets)

    def test_attribute_qoi_names_is_equal(self):
        self.assertEqual(self.qoi_config.qoi_names,
                         self.qoi_names)
