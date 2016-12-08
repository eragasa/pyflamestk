import pyflamestk.pyposmat
import pyflamestk.potentials

mc_sampler = pyflamestk.pyposmat.MonteCarloParameterSampler()
mc_sampler.create_lammps_simulations()
mc_sampler.run()
