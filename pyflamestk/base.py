import copy
import subprocess
import numpy as np
    
class Atom:
   def __init__(self, symbol, position):
       self._symbol = symbol
       self._position = position
       self._magnetic_moment = 0
       
   @property
   def symbol(self): 
       return self._symbol
       
   @symbol.setter
   def symbol(self,s): 
       self._symbol = s
       
   @property
   def position(self): 
       return np.array(self._position)
       
   @position.setter
   def position(self,p):
       self._position = np.aray(p)
           
   @property
   def magnetic_moment(self): 
       return self._magnetic_moment
       
   @magnetic_moment.setter
   def magnetic_moment(self,m): 
       self._magnetic_moment = m

class Structure:
    
    def __init__(self, obj=None):
        if obj is None:
            self._noncopy_init()
        else:
            self._copy_init(obj)
        self._ptol = 1.e-5
        self._vacancies = []
        self._interstitials = []
            
    def _noncopy_init(self):
        self._atoms = []
        self._structure_comment = None
        self._lattice_parameter = 1.0
        self._h_matrix = np.zeros(shape=[3,3])
        
    def _copy_init(self, obj):
        self._structure_comment = obj.structure_comment
        self._lattice_parameter = obj.lattice_parameter
        self._h_matrix = np.array(obj.h_matrix)
        
        self._atoms = copy.deepcopy(obj.atoms)

    @property
    def a1(self):
        a0 = self.lattice_parameter
        h1 = self.h1
        a1 = a0 * h1.dot(h1)**0.5
        return a1

    @property
    def a2(self):
        a0 = self.lattice_parameter
        h2 = self.h2
        a2 = a0 * h2.dot(h2)**0.5
        return a2

    @property
    def a3(self):
        a0 = self.lattice_parameter
        h3 = self.h3
        a3 = a0*h3.dot(h3)**0.5
        return a3

    @property
    def lattice_parameter(self): 
        return self._lattice_parameter
       
    @lattice_parameter.setter
    def lattice_parameter(self, a):
        if a > 0:
            self._lattice_parameter = a
        else:
            raise ValueError

    @property
    def structure_comment(self): 
        return self._structure_comment
       
    @structure_comment.setter
    def structure_comment(self, s):
        self._structure_comment = s

    @property
    def h_matrix(self):
        return self._h_matrix

    @h_matrix.setter
    def h_matrix(self,h):
        self._h_matrix = np.array(h)

    @property
    def h1(self):
        return np.array(self._h_matrix[0,:])

    @h1.setter
    def h1(self, h1):
        self._h_matrix[0,:] = np.array(h1)

    @property
    def h2(self):
        return np.array(self._h_matrix[1,:])

    @h2.setter
    def h2(self,h2):
        self._h_matrix[1,:] = np.array(h2)
 
    @property
    def h3(self):
        return np.array(self._h_matrix[2,:])

    @h3.setter
    def h3(self,h3):
        self._h_matrix[2,:] = np.array(h3)

    @property
    def b1(self):
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)

        b1 = 2 * np.pi * np.cross(a2,a3) / np.dot(a1, np.cross(a2,a3))
        return b1

    @property
    def b2(self):
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)


        b2 = 2 * np.pi * np.cross(a3,a1) / np.dot(a2, np.cross(a3,a1))
        return b2

    @property
    def b3(self):
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)

        b3 = 2 * np.pi * np.cross(a1,a2) / np.dot(a3, np.cross(a1,a2))
        return b3
        
    @property
    def n_atoms(self):
        return(len(self._atoms))

    @property
    def symbols(self):
        symbols = []
        for a in self._atoms:
            if not(a.symbol in symbols):
                symbols.append(a.symbol)
        return list(symbols)

    @property
    def atoms(self):
        return self._atoms
        
    @atoms.setter
    def atoms(self, atoms):
        self._atoms  = copy.deepcopy(atoms)

    @property
    def vacancies(self):
        return self._vacancies

    @property
    def interstitials(self):
        return self._interstitials

    def check_if_atom_exists_at_position(self,symbol,position):
        is_atom_exists = False           # initialize return variable

        # check to see if atom exists
        for a in self._atoms:
            diff = [abs(position[i]-a.position[i]) for i in range(3)]
            if max(diff) < self._ptol:
                print('new:  ',symbol,position)
                print('exist:',a.symbol,a.position)
                is_atom_exists = True
                return is_atom_exists # = True
                
        return is_atom_exists # = False

    def add_atom(self, symbol, position):
        is_atom_exists = self.check_if_atom_exists_at_position(symbol,position)
        if is_atom_exists is True:
            err_msg = "Tried to add {} @ {} an atom already there"
            err_msg = err_msg.format(symbol,position)
            raise ValueError(err_msg)
        else:
            self._atoms.append(Atom(symbol,position))
            self._interstitials.append([symbol,position])

    def remove_atom(self,symbol, position):
        for i,a in enumerate(self._atoms):
            if (a.symbol == symbol):
                diff = [abs(position[j]-a.position[j]) for j in range(3)]
                if max(diff) < self._ptol:
                    self._atoms.remove(self._atoms[i])
                    self._vacancies.append([symbol,position])
                    return
        # no atom found
        err_msg = "Tried to remove {} @ {}, no atom found"
        err_msg = err_msg.format(symbol,position)
        raise ValueError(err_msg)

    def get_number_of_atoms(self, symbol=None):
        if symbol is None:
            n_atoms = len(self._atoms)
        else:
            n_atoms = 0
            for atom in self._atoms:
                if (atom.symbol == symbol):
                    n_atoms += 1
        return n_atoms

    def normalize_h_matrix(self):
        # change the h_matrix where the lattice parameter is 1.0
        for i in range(3):
            self._h_matrix[i,:] = self._h_matrix[i,:] * self._lattice_parameter
        self._lattice_parameter = 1.0

        # calculate the norm of the h1 vector
        norm_h1 = np.sqrt(self._h_matrix[i,:].dot(self._h_matrix[i,:]))
        # the norm of the h1 vector is the lattice parameter
        self._lattice_parameter = norm_h1

        # normalize the h-matrix by the norm of h1
        for i in range(3):
            self._h_matrix[i,:] = self._h_matrix[i,:] / norm_h1

    def set_lattice_parameter(self, a1):
        # change the h_matrix where the lattice parameter is 1.0
        for i in range(3):
            self._h_matrix[i,:] = self._h_matrix[i,:] * self._lattice_parameter


        self._lattice_parameter = a1
        for i in range(3):
            self._h_matrix[i,:] = self._h_matrix[i,:] / self._lattice_parameter

    def __str__(self):
        str_out = "a = {}\n".format(self._lattice_parameter)
        str_out += "atom list:\n"
        for a in self._atoms:
            str_t = "{} {:10.6f} {:10.6f} {:10.6f} {:10.6f}"
            str_out += str_t.format(a.symbol,
                                    a.position[0],
                                    a.position[1],
                                    a.position[2],
                                    a.magnetic_moment)
 
def make_super_cell(structure, sc):
    supercell = Structure()
    supercell.structure_comment = "{}x{}x{}".format(sc[0],sc[1],sc[2])

    # set lattice parameter
    supercell.lattice_parameter = structure.lattice_parameter   

    # set h_matrix
    h = np.zeros(shape=[3,3])
    for i in range(3):
        h[i,:] = structure.h_matrix[i,:] * sc[i]
    supercell.h_matrix = h

    # add supercell atoms
    for i in range(sc[0]):
        for j in range(sc[1]):
            for k in range(sc[2]):
                for atom in structure.atoms:
                    symbol = atom.symbol
                    position = atom.position
                    position = [(i+position[0])/sc[0],\
                                (j+position[1])/sc[1],\
                                (k+position[2])/sc[2]]
                    supercell.add_atom(symbol,position)

    # return a copy of the supercell
    return copy.deepcopy(supercell)

class Simulation:
    def __init__(self):
        raise NotImplementedError
    def run(self):
        raise NotImplementedError

# assumes a unix-like operating system
def tail(fname, n_lines, offset=None):
    cmd_str = "/usr/bin/tail -n {} {}".format(str(n_lines), fname)
    p = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    lines = stdout.decode('ascii').splitlines()
    return lines

