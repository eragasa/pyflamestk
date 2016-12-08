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
        return np.array(self._h_matrix)

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
        return copy.deepcopy(self._atoms)
        
    @atoms.setter
    def atoms(self, atoms):
        self._atoms = copy.deepcopy(atoms)

    def add_atom(self, symbol, position):
        self._atoms.append(Atom(symbol,position))

    def get_number_of_atoms(self, symbol=None):
        if symbol is None:
            n_atoms = len(self._atoms)
        else:
            n_atoms = 0
            for atom in self._atoms:
                if (atom.symbol == symbol):
                    n_atoms += 1
        return n_atoms
        
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
    supercell = Structure(structure)
    supercell.structure_comment = "{}x{}x{}".format(sc[0],sc[1],sc[2])
    supercell.lattice_parameter = 1.
   
    h = np.zeros(shape=[3,3])
    for i in range(3):
        h[i,:] = structure.h_matrix[i,:] * sc[i] * structure.lattice_parameter

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
    return copy.deepcopy(supercell)

# assumes a unix-like operating system
def tail(fname, n_lines, offset=None):
    cmd_str = "/usr/bin/tail -n {} {}".format(str(n_lines), fname)
    p = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    lines = stdout.decode('ascii').splitlines()
    return lines

