#!/bin/env python
import unittest
import pyflamestk.potentials

class PotentialMorseWithSingleAtom(unittest.TestCase):

    def test_potentials_morse_correct_pairname_single_atom(self):
        name  = 'pair'
        atoms = ['Ni']
        pairname = "{}{}".format(atoms[0],atoms[0])

        morse = pyflamestk.potentials.Morse(name,atoms)
        result   = morse.pair

        print(pairname)
        print(result)

        self.assertEqual(pairname, result)

    def test_potentials_morse_correct_name_single_atom(self):
        name  = 'pair'
        atoms = ['Ni']        
        pairname = "{}{}".format(atoms[0],atoms[0])
        potname = "{}{}".format(name,pairname)

        morse = pyflamestk.potentials.Morse(name,atoms)
        result = morse.name
        self.assertEqual(potname, result)        
#class PotentialMorseWithTwoAtoms(unit.TestCase):
    
if __name__ == "__main__":
   unittest.main()