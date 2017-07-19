# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import pyflamestk.base as base
import numpy as np

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


class CutoffFunction(object):
    def __init__(self):
        pass

def cutoff_shifted_force(r,v,rcut):
    """

        r - numpy.array of distances
        v - numpy.array of energies
    Returns:
        sv - (numpy.array) force shifted potential
    """
    return sv

def cutoff_shifted_energy(r,v,rcut):
    """

        r - numpy.array of distances
        v - numpy.array of energies
        rc - (float) cutoff distance
    Returns:
        sv - (numpy.array) energy shifted potential
    """
    return sv

class ShiftedForceCutoff(CutoffFunction):
    def eval(self, r):
        pass

def func_cutoff(r,rcut,h):
    x = (r - rcut)/h
    phi = (x**4)/(1+x**4)
    # get index values where r > rcut
    r_gt_rcut = np.where(r >= rcut)
    phi[r_gt_rcut] = np.zeros(len(r_gt_rcut))
    return phi



class EamPotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'eam'
        self.param_dict = None
        self._param_names = None

        self._pair_function = None
        self._density_function = None
        self._embedding_function = None

    @property
    def pair_function(self):
        return self._pair_function

    @pair_function.setter
    def pair_function(self,func):
        self._pair_function = func

    @property
    def density_function(self):
        return self._density_function

    @density_function.setter
    def density_function(self, func):
        self._density_function = func

    @property
    def embedding_function(self):
        return self._embedding_function

    @embedding_function.setter
    def embedding_function(self, func):
        self._embedding_function = func

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names  = ['p.{}'.format(v) for v in self.pair_function.parameter_names]
        self._param_names += ['d.{}'.format(v) for v in self.density_function.parameter_names]
        self._param_names += ['e.{}'.format(v) for v in self.embedding_function.parameter_names]

    #def make_setfl_file(self):
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
        self._param_names = []
        for i in range(n_symbols):
            for j in range(n_symbols):
                if i <= j:
                    s_i = symbols[i]
                    s_j = symbols[j]
                    self._param_names.append("{}{}_A".format(s_i,s_j))
                    self._param_names.append("{}{}_rho".format(s_i,s_j))
                    self._param_names.append("{}{}_C".format(s_i,s_j))

class EamElectronDensityFunction(Potential):
    pass

class ExponentialDensityFunction(EamElectronDensityFunction):
    def __init__(self,symbols):
        EamElectronDensityFunction.__init__(self,symbols)
        self._pot_type = 'elec_dens_exp'
        self.param_dict = None
        self._param_names

        self._determine_parameter_names()
        self._create_param_dictionary()

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self._determine_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names = []
        for s in self._symbols:
            self._param_names.append('{}_rho0'.format(s))
            self._param_names.append('{}_beta'.format(s))
            self._param_names.append('{}_r0'.format(s))

    def evaluate(self,r,symbol,params,rcut=0,h=1):
        err_msg = "cannot find {} parameter for {}"

        rho0 = None
        beta = None
        r0 = None

        if '{}_rho0'.format(symbol) in self._param_names:
            rho0 = params['{}_rho0'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('rho0',symbol))
            raise ValueError(str_out)

        if '{}_beta'.format(symbol) in self._param_names:
            beta = params['{}_beta'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('beta',symbol))
            raise ValueError(str_out)

        
        if '{}_r0'.format(symbol) in self._param_names:
            r0 = params['{}_r0'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('r0',symbol))
            raise ValueError(str_out)

        val = rho0 * np.exp(-beta*(r/r0-1))

        if rcut != 0:
            val = val * func_cutoff(r,rcut,h)

        return val
class MorsePotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'pair_morse'
        self.param_dict = None
        self._param_names = None

        self._determine_parameter_names()
        self._create_param_dictionary()

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self.param_dict = {}
        for v in self._param_names:
            self.param_dict[v] = None

    def evaluate(self,r,pair,params, rcut=0, h=1.):
        """
        Args:

            rcut(float) - the cutoff function.  If set to 0, then no cutoff
        function is applied.
        """
        err_msg = "MorsePotential cannot find {} parameter for {},{} pair"

        # free_params = De, a, re
        D0 = None
        a = None
        r0 = None

        if '{}{}_D0'.format(pair[0],pair[1]) in self._param_names:
            D0 = params['{}{}_D0'.format(pair[0],pair[1])]
        elif '{}{}_D0'.format(pair[1],pair[0]) in self._param_names:
            D0 = params['{}{}_D0'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('D0',pair[0],pair[1])
            raise ValueError(str_out)

        if '{}{}_a'.format(pair[0],pair[1]) in self._param_names:
            a = params['{}{}_a'.format(pair[0],pair[1])]
        elif '{}{}_a'.format(pair[1],pair[0]) in self._param_names:
            a = params['{}{}_a'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('a',pair[0],pair[1])
            raise ValueError(str_out)

        if '{}{}_r0'.format(pair[0],pair[1]) in self._param_names:
            r0 = params['{}{}_r0'.format(pair[0],pair[1])]
        elif '{}{}_r0'.format(pair[1],pair[0]) in self._param_names:
            r0 = params['{}{}_r0'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('r0',pair[0],pair[1])
            raise ValueError(str_out)

        val = D0 * ((1 - np.exp(-a*(r-r0)))**2 -1)
        #val = D0 * (1 - np.exp(-a*(r-r0)))**2

        if rcut != 0:
            print('using rcut')
            val = val * func_cutoff(r,rcut,h)

        return val

    def _determine_parameter_names(self):
        symbols = self._symbols
        n_symbols = len(symbols)
        self._param_names = []
        for i in range(n_symbols):
            for j in range(n_symbols):
                if i <=j:
                    s_i = symbols[i]
                    s_j = symbols[j]
                    self._param_names.append("{}{}_D0".format(s_i,s_j))
                    self._param_names.append("{}{}_a".format(s_i,s_j))
                    self._param_names.append("{}{}_r0".format(s_i,s_j))

class EamEmbeddingFunction(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = None
        self.param_dict = None
        self._symbols = [s for s in symbols]

class UniversalEmbeddingFunction(EamEmbeddingFunction):
    def __init__(self,symbols):
        EamEmbeddingFunction.__init__(self,symbols)
        self._pot_type = 'embed_universal'
        self._determine_parameter_names()
        self._create_param_dictionary()

    @property 
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self.param_dict = {}
        for v in self._param_names:
            self.param_dict[v] = None


    def _determine_parameter_names(self):
        self._param_names = []
        symbols = self._symbols
        for s in symbols:
            self._param_names.append("{}_F0".format(s))
            self._param_names.append("{}_p".format(s))
            self._param_names.append("{}_q".format(s))
            self._param_names.append("{}_F1".format(s))

    def eval(rho, symbol, params):

        F0 = None
        p = None
        q = None
        F1 = None

        if "{}_F0".format(symbol) in self._param_names:
            F0 = params["{}_F0".format(symbol)]
        if "{}_p".format(symbol) in self._param_names:
            p = params["{}_p".format(symbol)]
        if "{}_q".format(symbol) in self._param_names:
            q = params["{}_q".format(symbol)]
        if "{}_F1".format(symbol) in self_param_names:
            F1 = params["{}_F1".format(symbol)]

        val = F0 * ( q/(q-p)*rho**p - p/(q-p)*rho**q) + F1 *rho
        return val

class BjsEmbeddingFunction(EamEmbeddingFunction):
    def __init__(self,symbols):
        EamEmbeddingFunction.__init__(self.symbols)
        self._pot_type = 'bjs'

    @property
    def parameter_names(self):
        self._determin_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names = None
        sybols = self._symbols
        for s in range(symbols):
            self._param_names.append("{}_F0").format(s)
            self._param_names.append("{}_gamma").format(s)
            self._param_names.append("{}_F1").format(s)

    def eval(rho, symbol, params):

        F0 = None # free paramter
        gamma = None # free parameter
        F1 = None # free parameter

        if "{}_F0".format(symbol) in self._param_names:
            F0 = params["{}_F0".format(symbol)]
        if "{}_gamma".format(symbol) in self._param_names:
            p = params["{}_gamma".format(symbol)]
        if "{}_F1".format(symbol) in self_param_names:
            F1 = params["{}_F1".format(symbol)]

        val = F0*(1-gamma*np.ln(rho))*rho**gamma + F1*gamma
        return val

