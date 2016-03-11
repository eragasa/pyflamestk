import pyflamestk.pyposmat
import pyflamestk.potentials

mc_sampler = pyflamestk.pyposmat.MonteCarloParameterSampler()
mc_sampler.run()

#for key in pot_config.elements:
#    print(key,'charge:',pot_config.elements[key]['charge'])
#
#for key in pot_config.pair:
#    print(key, pot_config.pair[key]['type'])
#
#param_list = pot_config.get_param_list()
#print(param_list)
#for param in param_list:
#    print(param,":",pot_config.get_param_value(param))
#
#    
#import numpy as np
#n_simulations = 1000
#param_set = []
#for i_simulation in range(n_simulations):
#    new_param_set = np.zeros(shape=[1,len(param_list)])
#    
#    i = 0
#    for param in param_list:
#        a = pot_config.get_param_value(param)[1]
#        b = pot_config.get_param_value(param)[2]
#        new_param_set[0,i] = np.random.uniform(low = a, high = b)
#        i += 1
#    param_set.append(new_param_set[0,:].tolist())
#np_param_set = np.array(param_set)
#print(np_param_set.shape)


#buckingham = pyflamestk.potentials.Buckingham(name = name,
#                                              atoms = atoms)