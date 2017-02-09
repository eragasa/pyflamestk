import pyflamestk.base as base
import pyflamestk.vasp as vasp
import numpy as np

class RocksSubmissionScript:
    def __init__(self,
                 fname,
                 job_name,
                 n_processors = 24,
                 queue_names = []):
        self._fname = fname

    def write(self,fname = None):
        if fname != None:
            self._fname = fname


class Simulation:
    def __init__(self, job_name, dir_name = None):
        self._job_name      = job_name
        
        # simulation directory
        if dir_name != None:
            self._dir_name = dir_name
        else:
            self._dir_name = self._job_name
        # input files
            
        self._fname_incar   = "INCAR"
        self._fname_poscar  = "POSCAR"
        self._fname_potcar  = "POTCAR"
        self._fname_contcar = "KPOINTS"

        # output files
        self._fname_outcar  = "OUTCAR"
        self._fname_oszicar = "OSZICAR"
        self._fname_wavecar = "WAVECAR"
        self._fname_contcar = "CONTCAR"

        # submission script
        self.fname_runjob = "runjob.job"

        self._is_done = False

    @property
    def is_done(self):
        base.tail(fname=self._fname_outcar)
        
    def _check_potcar_file(self):
        self._potcar = vasp.Potcar()
        self._potcar.read()
        print(self._potcar)

    def run(self):
        self._check_potcar_file()

    def is_simulation_complete(self):
        base.tail(self._fname_outcar)
        return True

class KpointsConvergenceTest:
    def __init__(self,
                 dir_potcar_lib,
                 fname_poscar="POSCAR",
                 fname_potcar=None,
                 fname_incar=None,
                 fname_log=None):

         self._lib_potcar = dir_potcar_lib
         self._fname_poscar = fname_poscar
         self._fname_potcar = fname_potcar
         self._fname_incar = fname_incar
         self._fname_log = fname_log
         self._outname_outcar = None
         self._logger = None
         
         self.poscar = None

         self._simulations = []

    def _log(self,string):
        print(string)
        if self._logger != None:
            self._logger.log(string)

    def _determine_kpoints_range(self):
        self._log("determining kpoints range from poscar file...")
        self._load_poscar_file()
        self.n_atoms = self.poscar.n_atoms
        self._log("n_atoms = {}".format(self.n_atoms))

        self.b1_norm = np.dot(self.poscar.structure.b1,
                              self.poscar.structure.b1) ** 0.5
        self.b2_norm = np.dot(self.poscar.structure.b2,
                              self.poscar.structure.b2) ** 0.5
        self.b3_norm = np.dot(self.poscar.structure.b3,
                              self.poscar.structure.b3) ** 0.5
    
        for l in range(3,10):
            self.k_point_mesh = [max(1,l + self.b1_norm + 0.5),
                                 max(1,l + self.b2_norm + 0.5),
                                 max(1,l + self.b2_norm + 0.5)]
            self._log("{} {}".format(l,self.k_point_mesh))
 
    def _load_poscar_file(self):
        self._log("poscar_file_name = {}".format(self._fname_poscar))
        self.poscar = vasp.Poscar(fname=self._fname_poscar)
        self.poscar.load(fname=self._fname_poscar)
        
    def _check_outcar_files(self):
        base.tail(self._fname_outcar)

    def run(self):
        self._log("running kpoints convergence test...")
        self._determine_kpoints_range()

kp_test = KpointsConvergenceTest('fuck_you')
kp_test.run()

sim = Simulation()
sim.run()