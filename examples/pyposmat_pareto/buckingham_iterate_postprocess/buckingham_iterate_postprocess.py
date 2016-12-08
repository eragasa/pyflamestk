#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
#from matplotlib.patches import cm
import numpy as np
import copy

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
 
pct_kept_each_iteration = 50           

#TODO:  do this programmatically
n_iterations = 10

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
        
                  
ps_aggregator = SimulationResultsAggregator(n_iterations=n_iterations)
ps_aggregator.report()

for p in ps_aggregator.param_list:
    if p != 'sim_id':
        p_mu_idx = ps_aggregator.param_evol_names.index("{}_mu".format(p))
        p_sigma_idx = ps_aggregator.param_evol_names.index("{}_sigma".format(p))
        x_data = ps_aggregator.param_evol_values[:,0]
        p_mu   = ps_aggregator.param_evol_values[:,p_mu_idx]
        p_std  = ps_aggregator.param_evol_values[:,p_sigma_idx]
        plt.plot(x_data,p_mu)
        plt.plot(x_data,p_mu-p_std)
        plt.plot(x_data,p_mu+p_std)
        plt.ylabel(p)
        plt.xlabel('iteration')
        plt.show()

#%%


#%%

for x_i, x_label in enumerate(ps_aggregator.qoi_err_keys):
    for y_i, y_label in enumerate(ps_aggregator.qoi_err_keys):
        if x_i < y_i:
            x_idx   = ps_aggregator.pareto_set_idx.index(x_label)
            y_idx   = ps_aggregator.pareto_set_idx.index(y_label)
            x_data  = ps_aggregator.n_iterations*[0] # intialize
            y_data  = ps_aggregator.n_iterations*[0] # initialize
            for i in range(ps_aggregator.n_iterations):
                x_data[i] = ps_aggregator.pareto_sets[i][:,x_idx]
                y_data[i] = ps_aggregator.pareto_sets[i][:,y_idx]
            
            plt.figure()
            color=iter(cm.rainbow(np.linspace(0,1,ps_aggregator.n_iterations)))
            plot_handles = ps_aggregator.n_iterations*[0]
            for i in range(ps_aggregator.n_iterations):
                c=next(color)
                plot_handles[i] = plt.scatter(x_data[i],
                                              y_data[i],
                                              label='{:n}'.format(i),
                                              color=c,s=2)
                pareto_front = pyflamestk.pareto.pareto_frontier_2d(ps_aggregator.pareto_sets[i][:,x_idx],
                                                                    ps_aggregator.pareto_sets[i][:,y_idx],
                                                                    maxX = False, maxY= False)
                plt.plot(pareto_front[0],pareto_front[1],color=c,linewidth=2)

            plot_handles.append(plt.scatter(ps_aggregator.aggregate_pareto_set[:,x_idx],
                                            ps_aggregator.aggregate_pareto_set[:,y_idx],
                                            label='all',color='k'))
                                            
            xmin = 0
            xmax = np.percentile(ps_aggregator.aggregate_pareto_set[:,x_idx],95)
            ymin = 0
            ymax = np.percentile(ps_aggregator.aggregate_pareto_set[:,y_idx],95)
            
            print(xmin,xmax,ymin,ymax)
            plt.axis([xmin,xmax,ymin,ymax])
            plt.legend(handles=plot_handles)
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()
            
#%%
pyflamestk.pareto.make_histograms_parameters_combined(sim_results[0])
#%% simulation results_i

pyflamestk.pareto.make_histograms_parameters_combined(sim_results)
#pyflamestk.pareto.make_histograms_qois_combined(sim_results,qois=qoi_keys)
#pyflamestk.pareto.make_pareto_set_jpdf_plot(sim_results=sim_results,
#                                            param_list = param_list)
#pyflamestk.pareto.make_cull_set_jpdf_plot(sim_results=sim_results,
#                                          param_list = param_list)
#pyflamestk.pareto.make_simulations_jpdf_plot(sim_results=sim_results,
#                                             param_list = param_list)

# for making comparative 2D Pareto plots
# pyflamestk.pareto.make_2d_pareto_plots(sim_results=sim_results)



#%% pareto plots with external data
#TODO: need to generalize this

def make_2d_pareto_plots_with_external_data(sim_results,
                                            sim_results_extern,
                                            sim_results_iterate,
                                            i_iterate,
                                            labels_extern_short = "",
                                            labels_extern_long = "",
                                            show_dominated=True,
                                            show_pareto=True,
                                            show_cull=True):
    n = len(sim_results.qoi_keys)
    for i in range(n):
        for j in range(n):
            if i < j and i != j:
                x_label = sim_results.qoi_keys[i]
                y_label = sim_results.qoi_keys[j]
                print("{} {}".format(x_label,y_label))
    
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                
                # apply performance requirements
                plt.figure()
                plt_handles = []
                if show_dominated:
                    plt_handles.append(plt.scatter(sim_results.np_all_sims[:,x_idx],
                                                   sim_results.np_all_sims[:,y_idx],
                                                   label='dominated',
                                                   marker='.',s = 1))
                if show_pareto:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set[:,x_idx],
                                                   sim_results.np_pareto_set[:,y_idx],
                                                   color='y',
                                                   label='Pareto'))
                if show_cull:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set_cull[:,x_idx],
                                                   sim_results.np_pareto_set_cull[:,y_idx],
                                                   color='r',
                                                   label='Pareto w/ constraints'))
                                                   
                if show_pareto:
                    pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                        sim_results.np_all_sims[:,y_idx],
                                                                        maxX = False, maxY= False)                
                    plt.plot(pareto_front[0],
                             pareto_front[1],
                             color='y',linewidth=2)
                             
                if show_cull:
                    pareto_front_cull = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_pareto_set_cull[:,x_idx],
                                                                             sim_results.np_pareto_set_cull[:,y_idx],
                                                                             maxX = False, maxY= False)
                    plt.plot(pareto_front_cull[0],
                             pareto_front_cull[1],
                             color='r',linewidth=2)
                             
                # show external
                x_idx_extern   = sim_results_extern.all_names.index(x_label)
                y_idx_extern   = sim_results_extern.all_names.index(y_label)
                

                plt.text(sim_results_extern.np_all_sims[0,x_idx_extern],
                         sim_results_extern.np_all_sims[0,y_idx_extern],
                         'LC2.0',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                plt.text(sim_results_extern.np_all_sims[1,x_idx_extern],
                         sim_results_extern.np_all_sims[1,y_idx_extern],
                         'BG2.0',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                plt.text(sim_results_extern.np_all_sims[2,x_idx_extern],
                         sim_results_extern.np_all_sims[2,y_idx_extern],
                         'BG1.7',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                #x_idx_iterate   = sim_results_iterate.all_names.index(x_label)
                #y_idx_iterate   = sim_results_iterate.all_names.index(y_label)

                #plt.text(sim_results_iterate.np_all_sims[i,x_idx_iterate],
                #         sim_results_iterate.np_all_sims[i,y_idx_iterate],
                #         "X",
                #         bbox={'facecolor':'green', 'alpha':1.0, 'pad':2}) 
                         
                            
                plt.axis([0,
                          1.5*max(np.percentile(sim_results.np_pareto_set[:,x_idx],90),
                                  max(sim_results_extern.np_all_sims[:,x_idx_extern])),
                          0,
                          1.5*max(np.percentile(sim_results.np_pareto_set[:,y_idx],90),
                                   max(sim_results_extern.np_all_sims[:,y_idx_extern]))])
                plt.legend(handles=plt_handles)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.show()         

# UNCOMMENT THIS OUT TO MAKE THIS WORK-------------------------------------------

fname_results_2 = "results_pub.out"          # original data                
sim_results_extern = pyflamestk.pareto.SimulationResults()
sim_results_extern.read_simulation_results(fname_in=fname_results_2)
sim_results_extern.set_qoi_type(qoi_type)

fname_results_3 = "iter_009/sim_results_009.dat"      # original data                
sim_results_iterate = pyflamestk.pareto.SimulationResults()
sim_results_iterate.read_simulation_results(fname_in=fname_results_3)
sim_results_iterate.set_qoi_type(qoi_type)


make_2d_pareto_plots_with_external_data(sim_results=ps_aggregator.sim_results[9],
                                        sim_results_extern=sim_results_extern,
                                        sim_results_iterate=sim_results_iterate,
                                        i_iterate = 2)

i_iter = 9
n_row, n_col = ps_aggregator.sim_results[9].np_pareto_set.shape


good_rows = []
for i_row in range(n_row):
    
    is_better = True
    
    str_out_err = ""
    for qoi in qoi_keys:
        qoi_err = "{}_{}".format(qoi,qoi_type)

        qoi_idx     = ps_aggregator.sim_results[i_iter].all_names.index(qoi_err)
        qoi_err_val = ps_aggregator.sim_results[i_iter].np_pareto_set[i_row,qoi_idx]

        qoi_idx_extern     = sim_results_extern.all_names.index(qoi_err)
        qoi_err_val_extern = sim_results_extern.np_all_sims[:,qoi_idx_extern].min()

        if qoi_err_val > qoi_err_val_extern:
            str_out_err += "\t{}: {:.4f} > {:.4f}\n".format(qoi_err,qoi_err_val,qoi_err_val_extern)
            is_better = False
        else:
            str_out_err += "\t{}: {:.4f} < {:.4f}\n".format(qoi_err,qoi_err_val,qoi_err_val_extern)
    print(i_row)
    print(str_out_err)
    
    str_out = "{:02} ".format(i_row)
    if is_better == True:
        good_rows.append(i_row)
        for p in param_list:
            p_idx = ps_aggregator.sim_results[9].all_names.index(p)
            p_val = ps_aggregator.sim_results[9].np_pareto_set[i_row,p_idx]

            #str_out += "{} ".format(p_val)
            
        for qoi in qoi_keys:
            qoi_idx = ps_aggregator.sim_results[i_iter].all_names.index(qoi)
            qoi_val = ps_aggregator.sim_results[i_iter].np_pareto_set[i_row,qoi_idx]
            str_out += "{:10.4f} ".format(qoi_val)
    else:
        str_out += "{:10}\n".format('FAILURE')
        str_out += str_out_err
    print(str_out)

#%%

# This code created plots for S. Phillpot's slides at MRS.
i_iter = 9
for ix, x in enumerate(qoi_keys):
    for iy, y in enumerate(qoi_keys):
        if x != y and ix < iy:

            # determine key for qoi error            
            x_key = "{}_{}".format(x,qoi_type)
            y_key = "{}_{}".format(y,qoi_type)
            
            # get x-axis pareto data
            x_i = ps_aggregator.sim_results[i_iter].all_names.index(x_key)        
            x_v = ps_aggregator.sim_results[i_iter].np_pareto_set[good_rows,x_i]

            # get y-axis pareto data
            y_i = ps_aggregator.sim_results[i_iter].all_names.index(y_key)        
            y_v = ps_aggregator.sim_results[i_iter].np_pareto_set[good_rows,y_i]

            # get comparative values
            x_idx_extern   = sim_results_extern.all_names.index(x_key)
            y_idx_extern   = sim_results_extern.all_names.index(y_key)

            plt_handles = []

            if qoi_type not in ['abserr','nabserr','sqerr','nsqerr']:
                plt_handles.append(plt.scatter(ps_aggregator.sim_results[9].qoi_values[x_key],
                                               ps_aggregator.sim_results[9].qoi_values[y_key],
                                               label='Experimental',
                                               color='k'))

            plt_handles.append(plt.scatter(x_v,y_v,color='y',label='Pareto'))                
            plt_handles.append(plt.scatter(sim_results_extern.np_all_sims[0,x_idx_extern],
                                           sim_results_extern.np_all_sims[0,y_idx_extern],
                                           label='LC+2.0',color='r'))
                       
            plt_handles.append(plt.scatter(sim_results_extern.np_all_sims[1,x_idx_extern],
                                           sim_results_extern.np_all_sims[1,y_idx_extern],
                                           label='BG+2.0',color='g'))

            plt_handles.append(plt.scatter(sim_results_extern.np_all_sims[2,x_idx_extern],
                                           sim_results_extern.np_all_sims[2,y_idx_extern],
                                           label='BG+1.7',color='b'))
            plt.xlabel(x_key)
            plt.ylabel(y_key)
            
            x_axis_min = min(x_v.min(),
                             sim_results_extern.np_all_sims[[0,1,2],x_idx_extern].min())
            x_axis_max = max(x_v.max(),
                             sim_results_extern.np_all_sims[[0,1,2],x_idx_extern].max())
            y_axis_min = min(y_v.min(),
                             sim_results_extern.np_all_sims[[0,1,2],y_idx_extern].min())
            y_axis_max = max(y_v.max(),
                             sim_results_extern.np_all_sims[[0,1,2],y_idx_extern].max())
            x_length = x_axis_max - x_axis_min
            y_length = y_axis_max - x_axis_min
            x_axis_max = x_axis_max + 0.1 * x_length
            y_axis_max = y_axis_max + 0.1 * y_length
            
            
            

            plt.legend(handles=plt_handles,loc='lower right')
            plt.axis([x_axis_min,x_axis_max,y_axis_min,y_axis_max])
            plt.show()
            
#%% 3D pareto plot

i_iter = 2
x_i = ps_aggregator.sim_results[i_iter].all_names.index('MgO_NaCl_G_abserr')        
x_v = ps_aggregator.sim_results[i_iter].np_pareto_set[good_rows,x_i]

# get y-axis pareto data
y_i = ps_aggregator.sim_results[i_iter].all_names.index('MgO_NaCl_B_abserr')        
y_v = ps_aggregator.sim_results[i_iter].np_pareto_set[good_rows,y_i]

# get y-axis pareto data
z_i = ps_aggregator.sim_results[i_iter].all_names.index('MgO_NaCl_latt_abserr')        
z_v = ps_aggregator.sim_results[i_iter].np_pareto_set[good_rows,z_i]

fig = plt.figure()
tri = Delaunay(np.array([x_v,y_v]).T)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x_v,y_v,z_v, triangles=tri.simplices, cmap=plt.cm.Spectral)
ax.scatter(x_v, y_v, z_v)
ax.set_xlabel('MgO_NaCl_G_abserr')
ax.set_ylabel('MgO_NaCl_B_abserr')
ax.set_zlabel('MgO_NaCL_latt')
ax.view_init(45, 45)
plt.show()

#--------------------- ---------------------------------------------------------
#%%

all_names = param_list + qoi_keys

all_names_idx = []
for name in all_names:
    all_names_idx.append(ps_aggregator.sim_results[i_iter].all_names.index(name))

good_pot = ps_aggregator.sim_results[i_iter].np_pareto_set[:,:]
good_pot = good_pot[good_rows,:]
good_pot = good_pot[:,all_names_idx]

n_rows, n_cols = good_pot.shape

print("n_rows: {}".format(n_rows))
print("n_cols: {}".format(n_cols))

#header string
str_out = ' '.join(all_names) + '\n'
i=0
for row in good_pot:
    row_str = ['{:10.6f}'.format(n) for n in row]
    row_str = ' '.join(row_str) + '\n'
    row_str = row_str.lstrip()
    str_out += row_str

print(str_out)

#%%
all_names_idx = []
for name in all_names:
    all_names_idx.append(sim_results_extern.all_names.index(name))

extern_pot = sim_results_extern.np_all_sims[:,:]
extern_pot = extern_pot[:,all_names_idx]
str_out = ' '.join(all_names) + '\n'
i=0
for row in extern_pot:
    row_str = ['{:10.6f}'.format(n) for n in row]
    row_str = ' '.join(row_str) + '\n'
    row_str = row_str.lstrip()
    str_out += row_str
    
print(str_out)
