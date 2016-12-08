#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
#from matplotlib.patches import cm
import numpy as np
import copy

class SimulationResultsAggregator:
    """
    This class aggregates several pareto sets.
    
    Atrributes:
    
    qoi_list (str[])
    param_list (str[])
    qoi_keys (str[])
    qoi_err_keys (str[])
    qoi_type (str)
    
    
    sim_results (obj array) - an python list of pyflamestk.pareto.SimulationResults()
        storing the data from each simulation result.
    pareto_set_idx (str[]) - a string of pareto set.
    pareto_sets (obj array) - an array of pareto sets colleccted from the 
        simulation results.  The objects are contained as numpy arrays with
        rows representing simulations.  Within the columns, values of parameters
        are followed by the simulation values of different quantities of
        interest, followed by the error type.  A list containing the list of
        idx values are contained within the attribute pareto_set_idx
    aggregate_pareto_set (numpy array) - the aggregated pareto set calculated
        from each of the individual pareto sets from each iteration.
    
    """
    def __init__(self,
                 param_list="",
                 qoi_list="",
                 qoi_type="abserr",
                 qoi_keys=None,
                 n_iterations=10):
        """
        
        Arguments:
        
        sim_results (pyflamestk.pareto.SimulationsResults[])
        param_list (str[]) - an array of strings, if left blank the parameter
            list will be taken from sim_results[0].param_list
        qoi_list (str[]) - an array of qoi of interests, if left blank the
            parameter list will be taken from sim_results[0].param_list
        """
        
        self.qoi_type = qoi_type
        self.qoi_err_type = qoi_type
        self.n_iterations = n_iterations
        self.__read_sim_results()
        # set qoi keys
        if qoi_keys == None:
            #TODO: This is a really stupid way of doing this because 
            #I really coded up pyflamestk.pareto.SimulationResults wrong
            #pyflamestk.pareto.SimulationResults should be fixed, but this
            #is a quick regex fix
        
            print("Setting the qoi keys from simulations results...")
            self.qoi_keys = copy.deepcopy(self.sim_results[0].qoi_keys)
            for i,s in enumerate(self.qoi_keys):
                new_key  = re.search(r'(.*)_(.*)',s).group(1)
                err_type = re.search(r'(.*)_(.*)',s).group(2)
                if err_type in ['err','abserr','nabserr','sqerr','nsqerr']:
                     self.qoi_keys[i] = new_key
        else:
            self.qoi_keys = qoi_keys
        
        # set the parameter list
        if param_list == "":
            self.param_list = self.sim_results[0].param_names.copy()
        else:
            self.param_list = param_list
            
        # set the qoi list
        if qoi_list == "":
            self.qoi_list = self.sim_results[0].qoi_names.copy()
        else:
            self.qoi_list = qoi_list
            
        self.pareto_sets = []              # intialize
        self.aggregate_dataset = []        # intialize
        self.aggregated_pareto_set = []    # intialize

        self.n_iterations = len(self.sim_results)            
        self.pareto_set_idx = self.param_list + self.qoi_list

        self.__aggregate_pareto_sets()
        self.__calculate_aggregate_pareto_set()
        self.__calculate_parameter_evolution()
        
    def report(self):
        self.get_number_of_datapoints()
        str_out = "Quantities of Interest:\n"            
        for k in self.qoi_keys:
            str_out += "\t {}\n".format(k)
        str_out +="Error Type: {}\n".format(self.qoi_err_type)
        str_out +="Number of Iterations Found; {}\n".format(self.n_iterations)
        str_out +="Number Of Simulations:\n"
        for i,n in enumerate(self.n_datapoints):
            str_out += "\t {} {}\n".format(i,n)
        str_out += "\t ttl {}\n".format(self.total_datapoints)
        print(str_out)

        print("Evolution of Parameter Values")
        for param in self.param_list:
            for i in range(self.n_iterations):
                p_idx   = self.sim_results[i].all_names.index(param)
                p_value = self.sim_results[i].np_all_sims[:,p_idx].mean()
                p_std   = self.sim_results[i].np_all_sims[:,p_idx].std()
                print("{} {} {} {}".format(i,param,p_value,p_std))
                
        print("Writing Evolution of Parameter Values")
        str_out = "iter_id"
        for p in self.param_list:
            if p != 'sim_id':
                str_out += "{}_mu {}_sigma ".format(p,p)
        str_out += "\n"
            
        for i in range(1, self.n_iterations):
            str_out += "{} ".format(i)
            for p in self.param_list:
                if p != 'sim_id':
                    p_idx   = self.sim_results[i].all_names.index(p)
                    p_value = self.sim_results[i].np_all_sims[:,p_idx].mean()
                    p_std   = self.sim_results[i].np_all_sims[:,p_idx].std()
                    str_out += "{} {} ".format(p_value,p_std)
            str_out += "\n"
        #print(str_out)        
            
    def __calculate_parameter_evolution(self):
        #create numpy array for evolution of numpy values
        self.param_evol_names = ['iter_id']
        for p in self.param_list:
            if p != 'sim_id':
                self.param_evol_names.append("{}_mu".format(p))
                self.param_evol_names.append("{}_sigma".format(p))

        ps_evolution = []
        for i in range(1, self.n_iterations):
            row = []
            row.append(i)
            for p in self.param_list:
                if p != 'sim_id':
                    p_idx   = self.sim_results[i].all_names.index(p)
                    p_value = self.sim_results[i].np_all_sims[:,p_idx].mean()
                    p_std   = self.sim_results[i].np_all_sims[:,p_idx].std()
                    print(p,p_idx,p_value,p_std)
                    row.append(p_value)
                    row.append(p_std)
            ps_evolution.append(row)
        
        self.param_evol_values = np.array(ps_evolution)
        

    def __read_sim_results(self):
        self.sim_results = []
        for idx in range(self.n_iterations):
            fname_results = "iter_{:03d}/sim_results_{:03d}.dat".format(idx,idx)
            fname_pareto  = "iter_{:03d}/pareto_results_{:03d}.dat".format(idx,idx)
            fname_cull    = "iter_{:03d}/cull_results_{:03d}.dat".format(idx,idx)            
            
            self.sim_results.append(pyflamestk.pareto.SimulationResults())
            self.sim_results[idx].read_simulation_results(fname_in=fname_results)
            self.sim_results[idx].set_qoi_type(qoi_type)
            
            #TODO: read this from configuration file
            self.sim_results[idx].qoi_values = {}
            self.sim_results[idx].qoi_values['MgO_NaCl_latt'] = 4.1212
            self.sim_results[idx].qoi_values['MgO_NaCl_B']    = 226.
            self.sim_results[idx].qoi_values['MgO_NaCl_G']    = 92.
            
            #TODO: if filename does not exist
            self.sim_results[idx].calculate_pareto_set(qoi_keys=qoi_keys)
            #sim_results[idx].write(set_type='pareto',fname='pareto.dat')
            
            #TODO: if filename does not exist
            self.sim_results[idx].cull_by_percentile(pct_kept=pct_kept_each_iteration)
            #sim_results[idx].write(set_type='cull',fname='cull.dat')

            self.sim_results[idx].calculate_parameter_estimates(param_list=param_list)
            self.sim_results[idx].calculate_qoi_estimates(qoi_keys=qoi_keys)
 
    def get_number_of_datapoints(self):
        self.n_datapoints = []
        self.total_datapoints = 0
        for i in range(self.n_iterations):
            (n_rows, n_cols) = self.pareto_sets[i].shape
            self.n_datapoints.append(n_rows)
            self.total_datapoints += n_rows
        
    def __aggregate_pareto_sets(self):
        
        self.pareto_sets = n_iterations*[0]        
        for i in range(len(sim_results)):
            col_vals = len(self.pareto_set_idx) * ['']
            for j,p in enumerate(self.pareto_set_idx):
                idx = self.sim_results[i].all_names.index(p)
                col_vals[j] = self.sim_results[i].np_pareto_set[:,idx]
            self.pareto_sets[i] = np.array(col_vals).T
            
    def __calculate_aggregate_pareto_set(self):

        # aggregate all the pareto sets into a single pareto set        
        self.all_pareto_sets = np.concatenate(self.pareto_sets,axis=0)
        (n_pareto, n_cols) = self.all_pareto_sets.shape

        # set an index
        self.all_pareto_sets[:,0] = np.arange(n_pareto,dtype=np.int16)

        self.__create_dataset_for_pareto_analysis()
        pyflamestk.pareto.bruteforce_algo(self.pareto_dataset)
    
        # mark pareto set
        pareto_set_ids = []      # initialize
        pareto_set = []          # initialize
        for s in self.pareto_dataset:
            if s.paretoStatus == 1:
                pareto_set_ids.append(s.id)
                pareto_set.append(s.vec)
        pareto_set = -np.array(pareto_set)
         
        self.aggregate_pareto_set = self.all_pareto_sets[pareto_set_ids]
        
    def __create_dataset_for_pareto_analysis(self, 
                                             qoi_keys="",
                                             qoi_type="abserr"):
        """
        
        Arguments:
      
        qoi_keys (str[]) - A list of the quantities of interest we are interested 
            in.  If qoi_keys are not set, then the qois contained in the first 
            pareto set (sim_results[0]) will be used.
        qoi_type (str[]) - The type of qoi error type which will be used in this 
            simulation.
            
          
        Returns:

        pareto_dataset (float[][]) - a regular number python array provided
            as a list of lists.
            
        Attributes:
        
        self.pareto_dataset - a regular number python array provided as a list
            of lists to use in pyflamestk.pareto.bruteforce_algo()
        """

        # Pareto set for all qois
        if qoi_keys == "":
            pass
        else:
            self.qoi_keys = copy.deepcopy(qoi_keys)
            
        self.qoi_err_keys = copy.deepcopy(self.qoi_keys)
        for i,s in enumerate(self.qoi_keys):
            if self.qoi_type == 'err':
                # TODO: implement this
                self.qoi_err_keys[i] = "{}_err".format(self.qoi_keys[i])
                raise ValueError("qoi_type: err has not been implemented yet")
            elif self.qoi_type == 'abserr':
                self.qoi_err_keys[i] = "{}_abserr".format(self.qoi_keys[i])
            elif self.qoi_type == 'nabserr':
                self.qoi_err_keys[i] = "{}_nabserr".format(self.qoi_keys[i])
            elif self.qoi_type == 'sqerr':
                self.qoi_err_keys[i] = "{}_sqerr".format(self.qoi_keys[i])
            elif self.qoi_type == 'nsqerr':
                self.qoi_err_keys[i] = "{}_nsqerr".format(self.qoi_keys[i])
            else:
                raise ValueError("unknown_qoi_type")
                  
        # check to see if each key exists
        for qoi_err_key in self.qoi_err_keys:
            if not qoi_err_key in self.qoi_list:
                errmsg = "qoi_error_key ({}) does not exist"
                errmsg = errmsg.format(qoi_err_key)
                raise ValueError(errmsg)
              
        # make dataset 
        self.pareto_dataset = []
        for n, sim in enumerate(self.all_pareto_sets):
            self.pareto_dataset.append(pyflamestk.pareto.Datapoint(n))
            for k in self.qoi_err_keys:
                idx = self.pareto_set_idx.index(k)
                self.pareto_dataset[n].addNumber(-sim[idx])
         
        return copy.deepcopy(self.pareto_dataset)
        
#------------------------------------------------------------------------------
        
fname_results_1 = "results.out"      # data
fname_out       = "params.out"       # data

# --- parameter info ---
param_list = ['chrg_Mg',  'chrg_O',  
              'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho', 
              'p_OO_a',   'p_OO_c',   'p_OO_rho',
              'p_MgO_a',  'p_MgO_c',  'p_MgO_rho']
# chrg_Mg = -chrg_O, no MgMg interaction
param_list_not_free = ['chrg_O', 
                       'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho',
                       'p_MgO_c']

qoi_type     = "abserr"
qoi_keys     = ['MgO_NaCl_latt', 'MgO_NaCl_B', 'MgO_NaCl_G']

n_iterations = 10 
pct_kept_each_iteration = 50           
ps_aggregator = SimulationResultsAggregator(n_iterations=n_iterations)
ps_aggregator.report()

#TODO:  do this programmatically

# do pareto analysis on iteration set

def write_pareto_set_to_file(sim_results, fname_out):
    """
    
    Arguments:
    
    sim_results (pyflamestk.pareto.SimulationResults)
    """
    err_msg = "write_pareto_set_to_file() not implemented"
    raise NotImplementedError(err_msg)
    

#%%
import re

free_param_list = ['chrg_Mg', 'p_OO_a', 'p_OO_rho', 'p_OO_c', 'p_MgO_a',  'p_MgO_rho']
qoi_err_type = qoi_type

    
def get_parameter_manifold(sim_result, param_list):
    param_list     = ["sim_id"] + param_list
    param_idx_list = [sim_result.all_names.index(p) for p in param_list]
    print("parameter list:\n",param_list)
    param_manifold = sim_result.np_pareto_set[:,param_idx_list]
    return param_list, param_manifold
    
def get_qoi_err_manifold(sim_result, qoi_keys, qoi_err_type):
    qoi_err_keys = ["{}_{}".format(k,qoi_err_type) for k in qoi_keys]
    qoi_err_keys = ['sim_id'] + qoi_err_keys 
    print("qoi error keys:\n",qoi_err_keys)
    qoi_err_idx_list = [sim_result.all_names.index(q) for q in qoi_err_keys]
    qoi_err_manifold = sim_result.np_pareto_set[:,qoi_err_idx_list]
    return qoi_err_keys, qoi_err_manifold
    
p_keys, p_manifold = get_parameter_manifold(ps_aggregator.sim_results[9],
                                            free_param_list)
e_keys, e_manifold = get_qoi_err_manifold(ps_aggregator.sim_results[9],
                                          qoi_keys,
                                          qoi_err_type)
                                          
from time import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from sklearn import manifold
from sklearn.utils import check_random_state

# Next line to silence pyflakes.
Axes3D

# Variables for manifold learning.
n_neighbors = 10
n_samples = 1000

# Create our sphere.
random_state = check_random_state(0)
p = random_state.rand(n_samples) * (2 * np.pi - 0.55)
t = random_state.rand(n_samples) * np.pi

# Sever the poles from the sphere.
indices = ((t < (np.pi - (np.pi / 8))) & (t > ((np.pi / 8))))
colors = p[indices]
x, y, z = np.sin(t[indices]) * np.cos(p[indices]), \
    np.sin(t[indices]) * np.sin(p[indices]), \
    np.cos(t[indices])

# Plot our dataset.
plt.suptitle("Manifold Learning with %i points, %i neighbors"
             % (1000, n_neighbors), fontsize=14)

ax = fig.add_subplot(251, projection='3d')
ax.scatter(x, y, z, c=p[indices], cmap=plt.cm.rainbow)
try:
    # compatibility matplotlib < 1.0
    ax.view_init(40, -10)
except:
    pass

sphere_data = np.array([x, y, z]).T
fig = plt.figure(figsize=(15, 8))

# Perform Locally Linear Embedding Manifold learning
methods = ['standard', 'ltsa', 'hessian', 'modified']
labels = ['LLE', 'LTSA', 'Hessian LLE', 'Modified LLE']

for i, method in enumerate(methods):
    t0 = time()
    trans_data = manifold\
        .LocallyLinearEmbedding(n_neighbors, 2,
                                method=method).fit_transform(sphere_data).T
    t1 = time()
    print("%s: %.2g sec" % (methods[i], t1 - t0))

    ax = fig.add_subplot(252 + i)
    plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
    plt.title("%s (%.2g sec)" % (labels[i], t1 - t0))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')

# Perform Isomap Manifold learning.
t0 = time()
trans_data = manifold.Isomap(n_neighbors, n_components=2)\
    .fit_transform(sphere_data).T
t1 = time()
print("%s: %.2g sec" % ('ISO', t1 - t0))

ax = fig.add_subplot(257)
plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
plt.title("%s (%.2g sec)" % ('Isomap', t1 - t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

# Perform Multi-dimensional scaling.
t0 = time()
mds = manifold.MDS(2, max_iter=100, n_init=1)
trans_data = mds.fit_transform(sphere_data).T
t1 = time()
print("MDS: %.2g sec" % (t1 - t0))

ax = fig.add_subplot(258)
plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
plt.title("MDS (%.2g sec)" % (t1 - t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

# Perform Spectral Embedding.
t0 = time()
se = manifold.SpectralEmbedding(n_components=2,
                                n_neighbors=n_neighbors)
trans_data = se.fit_transform(sphere_data).T
t1 = time()
print("Spectral Embedding: %.2g sec" % (t1 - t0))

ax = fig.add_subplot(259)
plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
plt.title("Spectral Embedding (%.2g sec)" % (t1 - t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

# Perform t-distributed stochastic neighbor embedding.
t0 = time()
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
trans_data = tsne.fit_transform(sphere_data).T
t1 = time()
print("t-SNE: %.2g sec" % (t1 - t0))

ax = fig.add_subplot(2, 5, 10)
plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
plt.title("t-SNE (%.2g sec)" % (t1 - t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

plt.show()