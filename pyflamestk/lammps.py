import subprocess
import shutil
import time
import os
import pyflamestk.base as base

class Simulation:

  def __init__(self, sim_name, directory, structure):
    self.name = sim_name
    self.fname_dir  = directory
    self.fname_structure = structure
    self.script = None
    self.potential = None

  def create(self):
    # make directory
    if os.path.exists(self.dir):
      pass
    else:
      os.path.makedirs(self.dir)

  def run(self):
    pass

class InputFile:
  def __init__(self):
    self.units = "metal"
    self.dimension = 3
    self.boundary = "p p p"
    self.atom_style = "atomic"
    self.atom_modify = "map array"
    self.structure_fname = ""

    self.potential_type = "eam/alloy"
    self.eam_alloy_file = "eam.alloy"
  def initializeSimulations(self):
    str_out  = "# written by pyPosMat\n"
    str_out += "# ---- intialize simulation\n"
    str_out += "clear\n"
    str_out += "units {}\n".format(self.units)
    str_out += "dimension {}\n".format(self.dimension)
    str_out += "boundary {}\n".format(self.boundary)
    str_out += "atom_style {}\n".format(self.atom_style)
    str_out += "atom_modify {}\n".format(self.atom_modify)
    return str_out

  def createAtoms(self):
    if self.structure_fname == "":
      print("LAMMPS structure file not specifiedi in creating LAMMPS file")
      exit()

    str_out  = "# ---- create atoms\n"
    str_out += "read_data {}\n".format(self.structure_fname)

  def defineEmpiricalPotential(self):
    str_out  = "# ---- define interatomic potential"
    if self.structure_fname == "eam/alloy":
      str_out += "pair_style eam/alloy\n"
      str_out += "pair_coeff * * {} Ni\n".format(self.eam_alloy_file)
      str_out += "neighbor 2.0 bin\n"
      str_out += "neigh_modify delay 10 check yes\n"
    return str_out

  def defineSettings(self):
    str_out  = "compute eng all pe/atom\n"
    str_out += "compute eatoms all reduce sum c_eng\n"
    return str_out

  def runMinimization(self):
    str_out  = "reset_timestep 0\n"
    str_out += "fix 1 all box/relax iso 0.0 vmax 0.001\n"
    str_out += "thermo 10\n"
    str_out += "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n"
    str_out += "min_style cg\n"
    str_out += "minimize 1e-25 1e-25 5000 10000\n"
    return str_out

  def defineVariables(self):
    str_out  = "# ---- define variables"
    str_out += "variable natoms equal \"count(all)\"\n"
    str_out += "variable tot_energy equal \"c_eatoms\"\n"
    str_out += "variable length_x equal \"lx/4\"\n"
    str_out += "variable length_y equal \"ly/4\"\n"
    str_out += "variable length_z equal \"lz/4\"\n"
    str_out += "variable ecoh equal \"pe/atoms\"\n"
    return str_out

  def output(self):
    str_out  = "# ---- output\n"
    str_out += "print \"pyPosMat output section\"\n"
    str_out += "print \"tot_energy = ${tot_energy}\"\n"
    str_out += "print \"num_atoms = ${natoms}\"\n"
    str_out += "print \"latt_const_a = ${length_x}\"\n"
    str_out += "print \"latt_const_b = ${length_y}\"\n"
    str_out += "print \"latt_const_c = ${length_z}\"\n"
    str_out += "print \"ecoh = ${ecoh}\"\n"
    str_out += "print \"lammps_sim_done\"\n"
    return str_out
  
class Structure(base.Structure):
   def __init__(self,obj):
       base.Structure.__init__(self,obj)

   def write(self, filename_out, symbol_list, atom_style):
      print("writing structure to {}".format(filename_out))

      xlo   = 0.0
      xhi   = self.h_matrix[0,0]
      ylo   = 0.0
      yhi   = self.h_matrix[1,1]
      zlo   = 0.0
      zhi   = self.h_matrix[2,2]
      xy    = self.h_matrix[0][1]
      xz    = self.h_matrix[0][2]
      yz    = self.h_matrix[1][2]

      file = open(filename_out,'w')
      # TODO: want to add more information such as sc, encoded on this first line.
      file.write("# {}\n".format(self.symbols))
      file.write("\n")
      file.write("{} atoms\n".format(self.get_number_of_atoms()))
      file.write("{} atom types\n".format(len(self.symbols)))
      file.write("\n")
      file.write("{:10.4f} {:10.4f} xlo xhi\n".format(xlo, xhi))
      file.write("{:10.4f} {:10.4f} ylo yhi\n".format(ylo, yhi))
      file.write("{:10.4f} {:10.4f} zlo zhi\n".format(zlo, zhi))
      file.write("\n")
      file.write("{:10.4f} {:10.4f} {:10.4f} xy xz yz\n".format(xy,xz,yz))
      file.write("\n")
      file.write("Atoms\n")
      file.write("\n")

      atom_id = 1
      for i_symbol, symbol in enumerate(symbol_list):
        for i_atom, atom in enumerate(self.atoms):
          if (atom.symbol == symbol):
            if atom_style == 'atomic':
              file.write("{} {} {:10.4f} {:10.4f} {:10.4f}\n".format(atom_id, 
                                                                     i_symbol + 1, 
                                                                     self.h_matrix[0,0]*atom.position[0],\
                                                                     self.h_matrix[1,1]*atom.position[1],\
                                                                     self.h_matrix[2,2]*atom.position[2]))
              atom_id += 1
            elif atom_style == 'charge':
                #TODO: the style of this file is atom_id, atom_type, charge, x, y, z
                # a dummy charge of q = 1 is put here, but might neeed to be improved in the future.
                q = 1
                file.write("{} {} {:10.4f} {:10.4f} {:10.4f} {:10.4f}\n".format(atom_id, 
                                                                     i_symbol + 1,
                                                                     q,
                                                                     self.h_matrix[0,0]*atom.position[0],\
                                                                     self.h_matrix[1,1]*atom.position[1],\
                                                                     self.h_matrix[2,2]*atom.position[2]))
                atom_id += 1
      file.close()

def make_lammps_structure_file(structure, fname_out, sym_order):
    structure_file = Structure(structure)
    structure_file.write(fname_out, sym_order)
# various helper functions

def checkLammpsSimulationDone(fname_lmps_log):
  lines = base.tail(fname = fname_lmps_log, n_lines = 1)
  if 'lammps_sim_done' in  lines:
    return True
  else:
    return False

def checkAllLammpsSimulationsDone(sim_directories):
  simIsDone = [False for i in sim_directories]
  for idx, dir in enumerate(sim_directories):
    fname = "{}/out.dat".format(dir)
    simIsDone[idx] = checkLammpsSimulationDone(fname)
    #print("{} {}".format(fname, simIsDone[idx]))
  if False in simIsDone:
    return False
  else:
    return True

def checkLammpsSimulationForError(fname_lmps_log):
  lines = pyflamestk.base.tail(fname = fname_lmps_log, n_lines = 1)
  for line in lines:
    if 'Neighbor list overflow, boost neigh_modify one' in line:
      return True
  else:
    return False

def checkAllLammpsSimulationsForErrors(sim_directories):
  simHasError = False
  for dir in sim_directories:
    fname = "{}/out.dat".format(dir)
    simHasError = checkLammpsSimulationForError(fname)
    if (simHasError):
      return True
  return False

def createLammpsSimulation(sim_name,
                           fname_structure,
                           fname_sim_template,
                           fname_potential = None,
                           fname_eam = None):
  dir_name = sim_name

  if os.path.exists(dir_name):
    shutil.rmtree(dir_name)

  shutil.copytree(fname_sim_template,
                  "{}".format(dir_name))

  shutil.copyfile(fname_structure,
                  "{}/lammps.structure".format(dir_name))

  if not(fname_potential == None):
    pass

  if not(fname_eam == None):
    shutil.copyfile(fname_eam,
                   "{}/eam.alloy".format(dir_name))

  for dir in sim_directories:
    p = subprocess.Popen("./runsimulation.sh", shell=True, cwd=dir)

  # check if simulations are done
  sims_finished = False
  while not(sims_finished):
    sims_finished = checkAllLammpsSimulationsDone(sim_directories = sim_directories)
    time.sleep(0.05)

def getCohesiveEnergy(fname_lammps_out):
  ecoh=0
  f = open(fname_lammps_out,'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("ecoh ="):
      ecoh = float(line.split('=')[1].strip())
  f.close()
  return ecoh

def getElasticComponent(fname_lammps_out, component):
  retval = 0
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("{} =".format(component)):
      retval = line.split('=')[1].strip().split(' ')[0]
      retval = float(retval)
  f.close()
  return retval

def getLatticeParameter(fname_lammps_out):
  alat = 0
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("latt_const_a"):
      alat = float(line.split('=')[1].strip())
  f.close()
  return alat

def getPressure(fname_lammps_out, type = 'total'):
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  pressure_xx = 0
  pressure_yy = 0
  pressure_zz = 0
  pressure_total = 0
  for line in lines:
    if line.startswith("pressure_total"):
      pressure_total = float(line.split('=')[1].strip())
    elif line.startswith("pressure_xx"):
      pressure_xx = float(line.split('=')[1].strip())
    elif line.startswith("pressure_yy"):
      pressure_yy = float(line.split('=')[1].strip())
    elif line.startswith("pressure_zz"):
      pressure_zz = float(line.split('=')[1].strip())
    else:
      pass
  if   type == 'total':
    return pressure_total
  elif type == 'xx':
    return pressure_xx
  elif type == 'yy':
    return pressure_yy
  elif type == 'zz':
    return pressure_zz
  elif type == 'all':
    return [pressure_xx, pressure_yy, pressure_zz, pressure_total]
  else:
    print("pyflamestk.lammps -> do not understand pressure_type: {}".format(type))

