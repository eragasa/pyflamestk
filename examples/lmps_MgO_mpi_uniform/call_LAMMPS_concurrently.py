import time

from mpi4py import MPI
from operator import itemgetter
import numpy as np
import subprocess, sys, os, shutil
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials

"""
This script runs a parameter sampling for the Buckingham potential using a
uniform distribution for each parameter in the parameter set and running
a lammps simulation.  This script does a uniform sampling of the parameter
ter space and uses MPI for parallelization.  Concurrency by providing a 
segregated file space for the workers.

Required Directories:
  lmp_scripts_db  
  structure_db

Required Files:
  potential.mod
  pyposmat.config
  pyposmat.potential
  pyposmat.qoi

Output Files:
  pyposmat.log - log file
  results.out - results file

Modify for the system:
  1.  lmp_scripts_db/*/runsimulation.sh provides the location of the lammps
  binary.  Alternatively, $LAMMPS_BIN variable can be set.
  2.  the pyflamestk module needs to be added to the $PYTHONPATH.
      export PYTHONPATH=$(cd ..;pwd):$PYTHONPATH
"""

# --- this class will eventually go into a pyflamestk package ---

class MpiPyPosmatUniformSampler(pyposmat.PyPosmatEngine2):

    def __init__(self,
                 fname_config_pyposmat = "pyposmat.config",
                 fname_config_potential = "pyposmat.potential",
                 fname_config_qoi = "pyposmat.qoi",
                 is_read = True,
                 is_restart = False, 
                 random_seed = None):

        self._project_dir = os.getcwd()        # set to current working dir
        self._worker_dir_pre = "mpi_worker"
        self._get_mpi_info()

        # call the constructor
        pyposmat.PyPosmatEngine2.__init__(self,
                fname_config_pyposmat = fname_config_pyposmat,
                fname_config_potential = fname_config_potential,
                fname_config_qoi = fname_config_qoi,
                is_read = is_read, 
                is_restart = is_restart,
                random_seed = random_seed)

        self.sampler_type = 'uniform'
        self._create_working_directory()

    @property
    def project_directory(self):
        return self._project_dir

    @project_directory.setter
    def project_directory(self, val):
        assert type(val), str
        self._project_dir = val

    @property
    def worker_directory_prefix(self):
        return self._worker_dir_pre

    @worker_directory_prefix.setter
    def worker_directory_prefix(self, val):
        assert type(val), str
        self._worker_dir_pre = val

    @property
    def mpi_rank(self): return self._rank

    @property
    def mpi_size(self): return self._size

    @property
    def mpi_name(self): return self._name

    def _get_mpi_info(self):
        self._comm = MPI.COMM_WORLD
        self._rank = self._comm.Get_rank()
        self._size = self._comm.Get_size()
        self._name = MPI.Get_processor_name()
         
    def _log(self,msg):
        """ 
        this method logs a message to a log file.  this overrides the original
        method by prepending the rank to the log message so that you can track
        messages on different MPI.

        Arguments:
        msg - a string which is the log message
        """
        log_msg = "[rank_{}] - {}".format(self._rank, msg)
        print(log_msg)
        self._f_log.write(log_msg+'\n')

    def _create_working_directory(self):
        """
        this method creates a working directory for each of the working
        processes.  This prevents conflict between the different working
        classes.   Have each worker move to their own directory so function 
        evaluations (e.g., LAMMPS files written) don't conflict.
        """
        print("project directory:{}".format(self._project_dir))
        # define the working directory
        self._worker_dir = os.path.join(self._project_dir, self._worker_dir_pre + "_" + str(self._rank))

        # if the path exists delete it
        if os.path.exists(self._worker_dir):
            msg_out = "DELETING existing Working directory:\n"
            msg_out += "\t" + self._worker_dir
            print(msg_out)
            sys.stdout.flush()
            shutil.rmtree(self._worker_dir)
        os.makedirs(self._worker_dir)

        # TODO: temp code
        shutil.copytree('lmp_scripts_db',os.path.join(self._worker_dir,'lmp_scripts_db'))
        shutil.copytree('structure_db',os.path.join(self._worker_dir,'structure_db'))

        os.chdir(self._worker_dir)

    def sample_parameter_space(self,
                               n_simulations):
        """
        This method samples a parameter space.  It currently overrides the
        inherited behavior.

        Arguments:
        n_simulations - the number of simulations (float)
        fname_results = the filename of the output file (default:'results.out')

        Returns:
        nothing returned
        """

        fname_results = "results_rnk_{}.out".format(self._rank)
        f = open(fname_results, 'w')

        start_sim_id = 0          # initialized to normally start at 0
        write_header_files = True # initialized to normally writer headers
        self._results = None      # initialize results
        
        # TODO: add restart code

        # header line strings
        header_line = ", ".join(['sim_id'] + self._names)
        types_line = ", ".join(['sim_id'] + self._types)

        # write headers
        if write_header_files:
            f.write(header_line + "\n")
            f.write(types_line + "\n")

        # do simulations
        for sim_id in range(start_sim_id,n_simulations):
            param_dict = None

            # choose sampler type
            if self._sampler_type == 'uniform':
                param_dict = self.get_uniform_sample()
            else:
                raise PyPosmatError('unknown sampler type')

            # evaluate simulations
            if sim_id % self._size == self._rank:
                n,t,v = self.evaluate_parameter_set(param_dict)
                data = [sim_id] + v
                self._log("r:{},sid:{}".format(self._rank,sim_id))
                data_line = ", ".join(str(s) for s in data)
                f.write(data_line + '\n')

                if self._results is None:
                    self._results = [data]
                else:
                    self._results.append(data)
        self._results = np.array(self._results)
        f.close()

        self._log("rank {} completed sims".format(self._rank))
        # wait until all MPI workers are done
        self._comm.Barrier()
        
        # clean up self
        os.chdir(self._project_dir)
        shutil.copy(os.path.join(self._worker_dir,fname_results),
                    os.path.join(self._project_dir,fname_results))
        shutil.rmtree(self._worker_dir)

        self._log("rank {} cleaned up workspace".format(self._rank))
        # wait for all to be done.
        self._comm.Barrier()

        # rank 0 finishes clean up
        if self._rank == 0:
            self._log("doing final clean")
            self._results = []

            # aggregate results
            for i in range(0,self._size):
                fname_results_in = "results_rnk_{}.out".format(i)
                names, types, results = self._read_results_file(fname_results_in)
                for result in results:
                    self._results.append(result)

            # sort by sim_id, column 0
            self._results = sorted(self._results, key=itemgetter(0))

            # write results
            f_out = open("results.out",'w')
            f_out.write(header_line + "\n")
            f_out.write(types_line + "\n")
            for result in self._results:
                str_line = ", ".join(str(s) for s in result) + "\n"
                f_out.write(str_line)
            f_out.close()

            # remove rank results
            for i in range(0,self._size):
                fname_results_in = "results_rnk_{}.out".format(i)
                os.remove(fname_results_in)

        # wait for rank 0 to clean everything up
        self._comm.Barrier()

# ---- begin script ----
random_seed = 0
n_simulations = 1000
mpi_uniform_sampler = MpiPyPosmatUniformSampler(is_restart = False, random_seed = random_seed)
mpi_uniform_sampler.sample_parameter_space(n_simulations)
