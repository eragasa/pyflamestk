# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import pyflamestk.base as base

class Potential(object):
    def __init__(self,symbols):
        self._pot_type = None
        self._symbols = symbols
        self._param_names = None

    @property
    def symbols(self):
        return self._symbols

    @symbols.setter
    def symbols(self,symbols):
        # type checking
        assert type(symbols),list
        for s in symbols:
            assert type(symbols),str

        self._symbols = symbols

    @property
    def potential_type(self):
        return self._pot_type

    @property
    def parameter_names(self):
        raise NotImplementedError

    def write_lammps_potential_file(self):
        raise NotImplementedError


    def _get_mass(self,element):
        if element == 'Mg':
            return 24.305
        elif element == "O":
            return 15.999
        else:
            raise ValueError("element {} not in database".format(element))
          
    def _get_name(self,element):
        if element == "Mg":
            return 'magnesium'
        elif element == "O":
            return 'oxygen'
        else:
            raise ValueError('element {} not in database'.format(element))
 
class BuckinghamPotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'buckingham'
        self.param_dict = None
        self._param_names = None

        self._determine_parameter_names()
        self._create_param_dictionary()
    @property
    def parameter_names(self):
        """list of str: list of parameter names"""
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self.param_dict = {}
        for v in self._param_names:
            self.param_dict[v] = None

    def _determine_parameter_names(self):
        symbols = self._symbols
        self._param_names = ["chrg_{}".format(s) for s in symbols]
        n_symbols =  len(symbols)
        for i in range(n_symbols):
            for j in range(n_symbols):
                if i <= j:
                    s_i = symbols[i]
                    s_j = symbols[j]
                    self._param_names.append("{}{}_A".format(s_i,s_j))
                    self._param_names.append("{}{}_rho".format(s_i,s_j))
                    self._param_names.append("{}{}_C".format(s_i,s_j))



