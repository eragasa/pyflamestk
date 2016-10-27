#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

fname_results_1 = "results_10000b.out"
qoi_type      = "abs_err"
qoi_list = ['MgO_NaCl_latt_abserr',
            'MgO_NaCl_c11_abserr',
            'MgO_NaCl_c12_abserr',
            'MgO_NaCl_c44_abserr',
            'MgO_NaCl_B_abserr',
            'MgO_NaCl_G_abserr']
            
param_list = ['chrg_Mg',  'chrg_O',     'p_MgMg_a',
              'p_MgMg_c', 'p_MgMg_rho', 'p_OO_a',
              'p_OO_c',   'p_OO_rho',   'p_MgO_a',
              'p_MgO_c',  'p_MgO_rho']
              
# TODO: fix this magic number
n_resamples = 10000


sim_results = pyflamestk.pareto.SimulationResults()
sim_results.read_simulation_results(fname_in=fname_results_1)
sim_results.set_qoi_type(qoi_type)
sim_results.create_dataset_for_pareto_analysis()
sim_results.calculate_pareto_set()

sim_results.add_performance_constraint('MgO_NaCl_latt_abserr', 0.03)
sim_results.add_performance_constraint('MgO_NaCl_B_abserr', 10.)
sim_results.add_performance_constraint('MgO_NaCl_G_abserr', 10.)
sim_results.apply_performance_constraints()

#%% MAKE 2D PARETO PLOTS
#n = len(qois)

n = len(sim_results.qois)
for i in range(n):
    for j in range(n):
        if i < j and i != j:
            x_label = sim_results.qois[i]
            y_label = sim_results.qois[j]
            print("{} {}".format(x_label,y_label))
            x_idx   = sim_results.all_names.index(x_label)
            y_idx   = sim_results.all_names.index(y_label)
            
            # calculate the pareto front
            pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                              sim_results.np_all_sims[:,y_idx],
                                              maxX = False, maxY= False)
                                              
            # plot results
            plt.figure()
            plt.scatter(sim_results.np_all_sims[:,x_idx],
                        sim_results.np_all_sims[:,y_idx])
            plt.scatter(sim_results.pareto_set[:,sim_results.qois.index(x_label)],
                        sim_results.pareto_set[:,sim_results.qois.index(y_label)],color='y')
            plt.plot(pareto_front[0],pareto_front[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()

#%% MAKE HISTOGRAMS

def make_histograms(sim_results):
  pass
  
param_names_id = []
print("Parameter names:")
for idx, param_name in enumerate(sim_results.param_names):
  print("{} {}".format(idx, param_name))
  param_names_id.append(sim_results.all_names.index(param_name))
param_names_id.remove(param_names_id[0])
#param_names_id.remove(param_names_id[3])


np_pareto_set = sim_results.np_all_sims[sim_results.pareto_set_ids]
for name in sim_results.param_names:
  x_idx = sim_results.all_names.index(name)
  try:
      
    
    data = np_pareto_set[:,x_idx]
    density = scipy.stats.gaussian_kde(data)
    x = np.linspace(min(data), max(data), 1000)

    plt.figure()
    plt.subplot(211)
    plt.hist(np_pareto_set[:,x_idx],bins=40)
    plt.subplot(212)
    plt.plot(x, density(x))
    plt.xlabel(name)
    plt.show()
  except np.linalg.linalg.LinAlgError:
    print("LinAlgError in {}".format(name))


#%% calculate covariance matrix, calculate correlation matrix
    
cov_matrix = np.cov(np.transpose(np_pareto_set[:,param_names_id]))
cor_matrix = np.cov(np.transpose(np_pareto_set[:,param_names_id]))
for i in range(len(param_names_id)):
    for j in range(len(param_names_id)):
        cor_matrix[i,j] = cov_matrix[i,j]/(np.sqrt(cov_matrix[i,i]*cov_matrix[j,j]))
print(cor_matrix)

#%% make density plots

def make_density_plot(x_data,y_data,x_label,y_label):
#    H, xedges, yedges = np.histogram2d(x_data,y_data,bins=15)
  gridsize = 30
  plt.hexbin(x_data,y_data, gridsize=gridsize)
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  cb = plt.colorbar()
  cb.set_label("freq")
  plt.show()

def make_jpdf_plot(x_data,y_data,x_label,y_label):
  xmin = x_data.min()
  xmax = x_data.max()
  ymin = y_data.min()
  ymax = y_data.max()
  X, Y = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
  positions = np.vstack([X.ravel(), Y.ravel()])
  values = np.vstack([x_data,y_data])
  kernel = scipy.stats.gaussian_kde(values)
  Z = np.reshape(kernel(positions).T, X.shape)
  #plt.figure()
  plt.pcolor(X,Y,Z)
  #plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
  #          extent=[xmin, xmax, ymin, ymax], aspect='equal')
  plt.plot(x_data, y_data, 'k.', markersize=2)
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.axis([xmin, xmax,ymin,ymax])
  #plt.set_ylim([ymin, ymax])
  cb = plt.colorbar()
  cb.set_label("freq")
  plt.show()

#%%
             
for pn_i, x_label in enumerate(param_list):
  for pn_j, y_label in enumerate(param_list):
    if pn_i < pn_j and pn_i != pn_j:
      x_data  = np_pareto_set[:,sim_results.all_names.index(x_label)]
      y_data  = np_pareto_set[:,sim_results.all_names.index(y_label)]
      make_density_plot(x_data,y_data,x_label,y_label)
      try:
        make_jpdf_plot(x_data,y_data,x_label,y_label)
      except np.linalg.linalg.LinAlgError:
        print("Kernel Density Estimator produced LinAlgErr")


#%% Cull the pareto set by the performance requirements
performance_requirements = {}
performance_requirements['MgO_NaCl_latt_abserr'] = 0.05
performance_requirements['MgO_NaCl_B_abserr'] = 30
performance_requirements['MgO_NaCl_G_abserr'] = 50

n_sims, n_qoi = np_pareto_set.shape

#make a copy of the pareto set
np_pareto_set_cull = np.copy(np_pareto_set)
rows_to_delete = []

#cull the pareto set by the performance requirements
for idx in range(n_sims):
  ps = np_pareto_set_cull[idx,:]
  is_delete_row = False
  for qoi_name in performance_requirements.keys():
    qoi_idx = sim_results.all_names.index(qoi_name)
    if ps[qoi_idx] > performance_requirements[qoi_name]:
      is_delete_row = True 
  if is_delete_row:
    rows_to_delete.append(idx)

# delete rows not in the pareto set
np_pareto_set_cull = np.delete(np_pareto_set_cull,rows_to_delete,axis=0)
    
# comparative pareto plots
n = len(sim_results.qois)
for i in range(n):
    for j in range(n):
        if i < j and i != j:
            x_label = sim_results.qois[i]
            y_label = sim_results.qois[j]
            print("{} {}".format(x_label,y_label))

            x_idx   = sim_results.all_names.index(x_label)
            y_idx   = sim_results.all_names.index(y_label)
            
            pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                sim_results.np_all_sims[:,y_idx],
                                                                maxX = False, maxY= False)
                                                                
            pareto_front_cull = pyflamestk.pareto.pareto_frontier_2d(np_pareto_set_cull[:,x_idx],
                                                                     np_pareto_set_cull[:,y_idx],
                                                                     maxX = False, maxY= False)

            # apply performance requirements
            plt.figure()
            plt.scatter(sim_results.np_all_sims[:,x_idx],
                        sim_results.np_all_sims[:,y_idx])
            plt.scatter(sim_results.pareto_set[:,sim_results.qois.index(x_label)],
                        sim_results.pareto_set[:,sim_results.qois.index(y_label)],color='y')
            plt.scatter(np_pareto_set_cull[:,x_idx],
                        np_pareto_set_cull[:,y_idx],color='r')
            plt.plot(pareto_front[0],pareto_front[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()

            # apply performance requirements
            plt.figure()
            plt.scatter(np_pareto_set_cull[:,x_idx],
                        np_pareto_set_cull[:,y_idx],color='y')
            plt.plot(pareto_front_cull[0],pareto_front_cull[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()
#%%   
           
for pn_i, x_label in enumerate(param_list):
  for pn_j, y_label in enumerate(param_list):
    if pn_i < pn_j and pn_i != pn_j:
      x_data  = np_pareto_set_cull[:,sim_results.all_names.index(x_label)]
      y_data  = np_pareto_set_cull[:,sim_results.all_names.index(y_label)]
      make_density_plot(x_data,y_data,x_label,y_label)
      try:
        make_jpdf_plot(x_data,y_data,x_label,y_label)
      except np.linalg.linalg.LinAlgError:
        print("Kernel Density Estimator produced LinAlgErr")

#%%----------------------------------------------------------------------
n = len(sim_results.qois)
for i in range(n):
    for j in range(n):
        if i < j and i != j:
            x_label = sim_results.qois[i]
            y_label = sim_results.qois[j]
            print("{} {}".format(x_label,y_label))
            x_idx   = sim_results.all_names.index(x_label)
            y_idx   = sim_results.all_names.index(y_label)
            pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                sim_results.np_all_sims[:,y_idx],
                                                                maxX = False, maxY= False)
            plt.figure()
            plt.scatter(sim_results.np_all_sims[:,x_idx],
                        sim_results.np_all_sims[:,y_idx])
            plt.scatter(sim_results.pareto_set[:,sim_results.qois.index(x_label)],
                        sim_results.pareto_set[:,sim_results.qois.index(y_label)],color='y')
            plt.plot(pareto_front[0],pareto_front[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()

#------------------------------------------------------------------------------
#%% parameter list

print("***culled_parameters***")
param_list = ['chrg_Mg',  'chrg_O',     
              'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho', 
              'p_OO_a',   'p_OO_c',   'p_OO_rho',
              'p_MgO_a',  'p_MgO_c',  'p_MgO_rho']

param_list_not_free = ['chrg_O', 'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho']
# remove parameters which aren't free
free_param_list = param_list.copy()    # make local copy to preserve param_list
free_param_list.remove('chrg_O')       # chrg_Mg = -chrg_O
free_param_list.remove('p_MgMg_a')     # no MgMg interaction
free_param_list.remove('p_MgMg_c')     # no MgMg interaction
free_param_list.remove('p_MgMg_rho')   # no MgMg interaction

# select only the parameters we are interested in
param_col_indxs = []                   # store indx of parameters
for param_name in free_param_list:
  param_indx = sim_results.all_names.index(param_name)
  print('{} {}'.format(param_name, param_indx))
  param_col_indxs.append(param_indx)
np_pareto_cull_params = np_pareto_set_cull[:,param_col_indxs]

#select only the parameters and do a KDE analysis
kernel = scipy.stats.gaussian_kde(np_pareto_cull_params.transpose())

#resample from the kernel
param_sample = kernel.resample(size=n_resamples)

#make comparative plots
for pn_i, x_label in enumerate(free_param_list):
  for pn_j, y_label in enumerate(free_param_list):
    if pn_i < pn_j and pn_i != pn_j:
      print("{}:{},{}:{}".format(pn_i,x_label,pn_j,y_label))
      
      x_data  = np_pareto_set[:,sim_results.all_names.index(x_label)]
      y_data  = np_pareto_set[:,sim_results.all_names.index(y_label)]
      try:
        make_jpdf_plot(x_data,y_data,x_label,y_label)
      except np.linalg.linalg.LinAlgError:
        print("Kernel Density Estimator produced LinAlgErr")
        
      # culled data set
      x1 = np_pareto_set_cull[:,sim_results.all_names.index(x_label)]
      y1 = np_pareto_set_cull[:,sim_results.all_names.index(y_label)]
      try:
        make_jpdf_plot(x1,y1,x_label,y_label)
      except np.linalg.linalg.LinAlgError:
        print("Kernel Density Estimator produced LinAlgErr")      
      # resampled data set
      x2 = param_sample.T[:,free_param_list.index(x_label)]
      y2 = param_sample.T[:,free_param_list.index(y_label)]
      try:
        make_jpdf_plot(x2,y2,x_label,y_label)
      except np.linalg.linalg.LinAlgError:
        print("Kernel Density Estimator produced LinAlgErr")      

#%% comparative 2d pareto plots

n = len(sim_results.qois)
for i in range(n):
    for j in range(n):
        if i < j and i != j:
            x_label = sim_results.qois[i]
            y_label = sim_results.qois[j]
            print("{} {}".format(x_label,y_label))
            
            x_idx   = sim_results.all_names.index(x_label)
            y_idx   = sim_results.all_names.index(y_label)
            pareto_front_all = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                    sim_results.np_all_sims[:,y_idx],                                                    
            plt.figure()
            plt.scatter(sim_results.np_all_sims[:,x_idx],
                        sim_results.np_all_sims[:,y_idx])
            plt.scatter(sim_results.pareto_set[:,sim_results.qois.index(x_label)],
                        sim_results.pareto_set[:,sim_results.qois.index(y_label)],color='y')
            plt.plot(pareto_front[0],pareto_front[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()

#%% write out parameter file
# This section is specifically for buckingham.  It has not been generalized
# for all potential applications

# get param_list_out
param_list_out = []
for idx, row in enumerate(param_sample.T):
  row_idx = idx
  row_params = []
  #print(idx,row)
  for param in param_list:
    param_value = ""
    if param in free_param_list:
      param_value = row[free_param_list.index(param)]
      #print("{} in free_param_list, value = {}".format(param,param_value))
    else:
      # assignment of non-free parameters
      if param == 'chrg_O':
        param_value = -row[free_param_list.index('chrg_Mg')]
      elif param == 'p_MgMg_a':
        param_value = 0
      elif param == 'p_MgMg_c':
        param_value = 0
      elif param == 'p_MgMg_rho':
        param_value = 0.5
      else:
        pass
    row_params.append(param_value)
  param_list_out.append(row_params)


# print header
str_out = ""
for param in param_list:
  str_out += "{} ".format(param)

# print params
for idx,row in enumerate(param_list_out):
  str_out += "\n"
  str_out += "{} ".format(idx)
  for param_value in row:
    str_out += "{} ".format(param_value)

fname_out = 'param_out.dat'
f = open(fname_out, mode='w')
f.write(str_out)  
f.close()

#%% READ PARAMETER OUT FILE

class ParameterFileReader:
  def __init__(self, fname_param_in = 'params.in'):
    self.fname_params = fname_param_in
    f = open(self.fname_params)
    self.lines = f.readlines()
    f.close()

    n_lines = len(self.lines)
    self.param_names = self.lines[0].strip().split(' ')
    self.params = []
    for idx in range(1,n_lines):
      self.params.append([float(num) for num in self.lines[idx].split()])
     
param_reader = ParameterFileReader(fname_param_in=fname_out)

#------------------------------------------------------------------------------
#%% make density plots

#norm_pareto_set = pareto_set
#norm_pareto_set[:,0] = norm_pareto_set[:,0]/0.02
#norm_pareto_set[:,1] = norm_pareto_set[:,1]/27.2
#norm_pareto_set[:,2] = norm_pareto_set[:,2]/39.5
#norm_pareto_set[:,3] = norm_pareto_set[:,3]/7.2
#norm_pareto_set[:,4] = norm_pareto_set[:,4]
#norm_pareto_set[:,5] = norm_pareto_set[:,5]

# plot violin plot
plt.violinplot(sim_results.pareto_set[:,0], showmeans=False, showmedians=True)
plt.violinplot(sim_results.pareto_set[:,1], showmeans=False, showmedians=True)
plt.violinplot(sim_results.pareto_set[:,2], showmeans=False, showmedians=True)
plt.violinplot(sim_results.pareto_set[:,3], showmeans=False, showmedians=True)
plt.violinplot(sim_results.pareto_set[:,4], showmeans=False, showmedians=True)
plt.violinplot(sim_results.pareto_set[:,5], showmeans=False, showmedians=True)

min_error = np.zeros((1,6))
min_error[:,1] = 1
min_error[:,2] = 1
min_error[:,3] = 1