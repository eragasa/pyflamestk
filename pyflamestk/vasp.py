import re
import copy
import pyflamestk.base as base
import numpy as np

class Incar():
    def __init__(self, fname="INCAR"):
        self._fnamename = fname

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, fn):
        self._filename = fn

class Outcar():
    def __init__(self, fname="OUTCAR"):
        self._filename = fname
        self._total_energy = None
        self._phonon_eig_val = None # phonon eigenvalues
        self._phonon_eig_vec = None # phonon eigenvectors
                
    @property
    def total_energy(self): return self._total_energy
        
    @property
    def phonons(self): 
        return self._phonon_eig_val, self._phonon_eig_vec
        
    def read(self):
        with open(self._filename) as f:
            while line in f:
                if "free  energy   TOTEN" in line:
                    E = line.strip().split('=')[1].strip().split(' ')[0]
                    E = float(E)
                    self._total_energy = E
                    
    def __string__():
        print("total_energy[eV] = {}".format(self._ecoh_per_structure))
    def get_phonons(self):
        pass
    
    def get_time_of_calculation(self):
        pass
    def get_number_of_atoms(self):
        pass
    
    def get_energy(self, fname="OUTCAR"):
        self.ecoh_per_structure = None
        self.ecoh_per_atom = None

class Potcar():
    def __init__(self, symbols = None,
                       potcar_dir = None,
                       fname="POTCAR"):                   
        self._potcar_dir = potcar_dir
        self._fname = fname
        
        self._symbols   = symbols
        self._encut_min = []
        self._encut_max = []
        self._models    = []
        self._exch      = []
    
    @property    
    def filename(self): return self._fname
        
    @filename.setter
    def filename(self, fn): self._fname = fn
        
    @property
    def min_energy_cutoff(self): min(self._encut_min)
        
    @property
    def max_energy_cutoff(self): max(self._encut_max)
        
    def read(self, fname = "POTCAR"):
        # initialize arrays
        self._symbols   = []
        self._encut_min = []
        self._enmin_max = []
        self._models = []
        self._xc      = []
        
        with open(fname,'r') as f:
            for line in f:
                line = line.strip()
                if 'TITEL' in line:
                    symbol = line.split('=')[1].strip().split(' ')[1]
                    self._symbols.append(symbol)
                elif 'ENMIN' in line:
                    enmin = float(line.split(';')[0].split('=')[1].strip())
                    enmax = float(line.split(';')[0].split('=')[1].strip())
                    self._encut_min.append(enmin)
                    self._encut_max.append(enmax)
                elif 'LEXCH' in line:
                    xc = line.split('=')[1].strip()
                    self._xc.append(xc)
                elif "VRHFIN" in line:
                    m = line.split('=')[1].strip()
                    self._models.append(m)
                
    def __str__(self):
        header_row   = "symbol enmin enmax xc\n"
        format_row   = "{}({}) {:10.6f} {:10.6f} {}\n"
        
        n_atoms = len(self._symbols)
        str_out      = header_row
        for i in range(n_atoms):
            str_out += format_row.format(self._symbols[i],
                                         self._models[i],
                                         self._encut_min[i],
                                         self._encut_max[i],
                                         self._xc[i])
        return str_out
    def write(self,fname = None):
        if fname != None:
            self._fame = fname
            
class Kpoints():
    def __init__(self, 
                 fname="KPOINTS",
                 comment_line="Automatic Mesh",
                 n_kpoints = 0,
                 mesh_type = "Monkhorst-Pack",
                 mesh_size = [4,4,4],
                 mesh_shift = [0,0,0]):
                     
        self._fname = fname
        self._comment_line = comment_line
        self._n_kpoints = n_kpoints
        self._mesh_type = mesh_type
        self._mesh_size = mesh_size
        self._mesh_shift = mesh_shift
        
    def to_string(self):
        str_out = "{}\n".self._comment_line
        str_out = "{}\n".self._n_kpoints
        str_out = "{}\n".self._mesh_type

        str_out = "{:3d} {:3d} {:3d}\n".format(self._mesh_size[0],
                                               self._mesh_size[1],
                                               self._mesh_size[2])
        str_out = "{:3d} {:3d} {:3d}\n".format(self._mesh_shift[0],
                                               self._mesh_shift[1],
                                               self._mesh_shift[2])
        return str_out
        
    def write(self):
        f = open(self._fname, 'w')
        f.write(self.to_string)
        f.close()

class Poscar(base.Structure):
    
    def __init__(self,obj=None):
        base.Structure.__init__(self,obj)             
        self._fname_in = 'POSCAR'
        self._fname_out = 'POSCAR.out'

    @property
    def filename_in(self):
        return self._fname_in

    @filename_in.setter
    def filename_in(self,fname_in):
        self._fname_in = fname_in

    @property
    def filename_out(self):
        return self.filename_out
        
    @filename_out.setter
    def filename_out(self,fname_out):
        self._fname_out = fname_out
    
    def read_file(self, fname_in=None):
        if fname_in is not None:
            self._fname_in = fname_in
            
        f = open(self._fname_in, 'r')
        
        self.structure_comment = f.readline().strip()
        self.lattice_parameter = float(f.readline())

        # read h_matrix
        h_matrix = np.zeros(shape=[3,3])
        for i in range(3):
            h_row = f.readline().strip().split()
            h_row = np.array([float(val) for val in h_row])
            h_matrix[i,:] = h_row
        self.h_matrix = h_matrix

        # read symbols and atoms per symbol
        symbols        = f.readline().strip().split()
        n_symbols_list = [int(i) for i in f.readline().strip().split()]

        # read in atoms
        coordinate_style = f.readline()[0].upper()
        assert coordinate_style in ['D','C']
        
        # direct coordinate type
        if coordinate_style == "D":
            for i_sym, sym in enumerate(symbols):
                n_atoms = n_symbols_list[i_sym]
                for i_atom in range(n_atoms):
                    line = f.readline().strip().split()
                    position = [float(line[i]) for i in range(3)]
                    self.add_atom(sym,position)
        # indirect coordinate type
        elif coordinate_style == "C":
            for i_sym, sym in enumerate(symbols):
                n_atoms = n_symbols_list[i_sym]
                for i_atom in range(n_atoms):
                    line = f.readline().strip().split()
                    position = [float(line[i]) for i in range(3)]
                    self.add_atom(sym,position)
        else:
            err_str = "unsupported vasp coordinate_type: {}".format(coordinate_style)
            raise ValueError(err_str)

    def write_file(self, fname_out):
        self._fname_out = fname_out
        
        comment_string = "automatically generated by pyflamestk.vasp.Poscar()\n"
        str_poscar  = comment_string
        str_poscar += "{:10.6}\n".format(self.lattice_parameter)
        
        # string for h-matrix
        for i in range(3):
            h_row_template = "{:10.6f} {:10.6f} {:10.6f}\n"
            str_poscar += h_row_template.format(self.h_matrix[i,0],
                                                self.h_matrix[i,1],
                                                self.h_matrix[i,2])
                                       
        sym_list = self.symbols
        str_atomlist = ""
        str_atomnum  = ""
        for sym in sym_list:
            nAtoms = self.get_number_of_atoms(sym)
            str_atomlist   += " " + sym
            str_atomnum    += " " + str(nAtoms)
        str_atomlist   += "\n"
        str_atomnum    += "\n"

        str_poscar += str_atomlist
        str_poscar += str_atomnum
        str_poscar += "Direct\n"

        for symbol in self.symbols:
            for i, atom in enumerate(self.atoms):
                if symbol == atom.symbol:
                    pos_template = "{:10.6f} {:10.6f} {:10.6f}\n"
                    str_position = pos_template.format(atom.position[0],
                                                       atom.position[1],
                                                       atom.position[2])
                    str_poscar += str_position
 
        f = open(fname_out, 'w')
        f.write(str_poscar)
        f.close()
        
def make_super_cell(obj, scp):
    sc = base.make_super_cell(copy.deepcopy(obj), list(scp))
    return copy.deepcopy(Poscar(sc))
