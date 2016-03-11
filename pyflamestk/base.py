import copy
import subprocess
class Atom:
   def __init__(self, symbol, position):
       self.symbol = symbol
       self.position = position
       self.magnetic_moment = 0

class Structure:
   def __init__(self):
      self.atomList = []
      self.structureComment = ''
      self.lattice_parameter = 1.0
      self.h_matrix = [[0 for x in range(3)] for x in range(3)]

   def addAtom(self, symbol, position):
      self.atomList.append(Atom(symbol,position))

   def getH1(self):
      return(self.h_matrix[0])

   def getH2(self):
      return(self.h_matrix[1])

   def getH3(self):
      return(self.h_matrix[2])

   def getNumberOfAtoms(self, symbol='all'):
     if (symbol == 'all'):
       return(len(self.atomList))
     else:
       atomCount = 0
       for atom in self.atomList:
         if (atom.symbol == symbol):
           atomCount = atomCount + 1
       return atomCount

   def getSymbolList(self):
     symbol_list = []
     for atom in self.atomList:
        if not(atom.symbol in symbol_list):
           symbol_list.append(atom.symbol)
     return copy.deepcopy(symbol_list)

   def print(self):
     print("lattice_parameter: {}".format(self.lattice_parameter))
     print("atoms:")
     for atom in self.atomList:
       print(atom.symbol, atom.position, atom.magnetic_moment)
 
   def makeSupercell(self, sc_parameters):
       supercell = Structure()
       supercell.structureComment = self.structureComment
       supercell.lattice_parameter = 1.

       for i in range(3):
           for j in range(3):
               supercell.h_matrix[i][j] = sc_parameters[i] * self.h_matrix[i][j] * self.lattice_parameter


       for i in range(sc_parameters[0]):
           for j in range(sc_parameters[1]):
               for k in range(sc_parameters[2]):
                   for atom in self.atomList:
                       symbol = atom.symbol
                       position = atom.position
                       position = [(i+position[0])/sc_parameters[0],\
                                   (j+position[1])/sc_parameters[1],\
                                   (k+position[2])/sc_parameters[2]]
                       supercell.addAtom(symbol,position)
       return copy.deepcopy(supercell)

# assumes a unix-like operating system
def tail(fname, n_lines, offset=None):
  cmd_str = "/usr/bin/tail -n {} {}".format(str(n_lines), fname)
  p = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  stdout, stderr = p.communicate()
  lines = stdout.decode('ascii').splitlines()
  return lines

