import pyflamestk.base

class FccUnitCell(pyflamestk.base.Structure):
  def __init__(self, a0, atom_sym):
    pyflamestk.base.Structure.__init__(self)

    self.lattice_parameter = a0
    for i in range(3):
      self.h_matrix[i][i] = 1.

    positions = [ [0.0, 0.0, 0.0],
                  [0.5, 0.5, 0.0],
                  [0.5, 0.0, 0.5],
                  [0.0, 0.5, 0.5] ]

    for position in positions:
      self.addAtom(atom_sym, position)

class BccUnitCell(pyflamestk.base.Structure):
  def __init__(self, a0, atom_sym):
    pyflamestk.base.Structure.__init__(self)

    self.lattice_parameter = a0
    for i in range(3):
      self.h_matrix[i][i] = 1.

    positions = [ [ 0.0, 0.0, 0.0],
                  [ 0.5, 0.5, 0.5] ]

    for position in positions:
      self.addAtom(atom_sym, position)

class HcpUnitCell(pyflamestk.base.Structure):
  def __init__(self, a0, atom_sym):
    pyflamestk.base.Structure.__init__(self)

    self.lattice_parameter = a0

    self.h_matrix = [ [ 1.00, 0.00000000, 0.000000],
                      [-0.50, 0.86602531, 0.000000],
                      [ 0.00, 0.00000000, 1.648811] ]

    positions = [ [1/3, 2/3, 1/4],
                  [2/3, 1/3, 3/4] ]

    for position in positions:
      self.addAtom(atom_sym, position)
          
class NaClUnitCell(pyflamestk.base.Structure):
    def __init__(self, a0, cation_sym = 'Na', anion_sym = 'Cl'):
        pyflamestk.base.Structure.__init__(self)
        
        self.lattice_parameter = a0
        for i in range(3):
            self.h_matrix[i][i] = 1.
            
        cation_positions = [ [ 1/2, 1/2,  0.  ],
                             [ 1/2, 0.,   1/2 ],
                             [ 0.,  0.,   0., ],
                             [ 0.,  1/2,  1/2 ]]
                             
        anion_positions  = [ [ 0.,  1/2,  0.   ],
                             [1/2,  0.,   0.   ],
                             [1/2,  1/2,  1/2  ],
                             [0.,   0.,   1/2  ]]
                             
        for position in cation_positions:
            self.addAtom(cation_sym, position)
            
        for position in anion_positions:
            self.addAtom(anion_sym, position)
            
            