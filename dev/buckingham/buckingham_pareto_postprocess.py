#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
import numpy as np

def pareto_frontier_2d(Xs, Ys, maxX=True, maxY=True):
    '''
    Method to take two equally-sized lists and return just the elements which
    lie on the Pareto frontier, sorted into order.  Default behaviour is to 
    find the maximum for both X and Y, but the option is available to specify 
    maxX = False or maxY = False to find the minimum for either or both of the 
    parameters.
    
    original code: Jamie Bull,
    '''
    # Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    # Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
    # Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
                p_front.append(pair) # … and add them to the Pareto frontier
        else:
            if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
                 p_front.append(pair) # … and add them to the Pareto frontier
    # Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

def read_pareto_set(fname_in='pareto.in'):
  pass

def write_pareto_set(param_names,
                     qoi_names,
                     pareto_set, 
                     fname_out='pareto.out'):
  f = open(fname_out,'w')
  str_out = ""
  for name in param_names:
      str_out += "{} ".format(name)
  str_out += "| "
  for name in qoi_names:
      str_out += "{} ".format(name)
  str_out += "\n"
  f.write(str_out)
  for myset in pareto_set:
      str_out = ""
      for item in myset:
          str_out += "{} ".format(item)
      str_out += "\n"
  f.write(str_out)
  f.close()
    
def make_pareto_plot(xlabel,ylabel):
    x_idx   = all_names.index(x_label)
    y_idx   = all_names.index(y_label)
    pareto_front = pareto_frontier_2d(np_all_sims[:,x_idx],
                                      np_all_sims[:,y_idx],
                                      maxX = False, maxY= False)
    plt.figure()
    plt.scatter(np_all_sims[:,x_idx],
                np_all_sims[:,y_idx])
    plt.scatter(pareto_set[:,i],
                pareto_set[:,j],color='y')
    plt.plot(pareto_front[0],pareto_front[1])
    plt.axis([min(pareto_front[0]),
              max(pareto_front[0]),
              min(pareto_front[1]),
              max(pareto_front[1])])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

def create_histogram(name):
    plt.figure()
    x_idx = all_names.index(name)
    plt.hist(np_all_sims[:,x_idx])
    plt.show()


class SimulationResults:
    
  def __init__(self):      
    self.performance_requirements = {}
    self.fname_pareto_out = "pareto.out"

    # initialize variables
    self.all_sims    = []     # python array of all simulation data
    self.param_names = []     # array of parameter names
    self.qoi_names   = []     # array of qoi names
    self.pareto_set  = []
    self.np_pareto_set = []
    self.np_all_sims = []     # numpy array of all simulation data
       
  def read_simulation_results(self,fname_in):

    f = open(fname_results_1)
    lines = f.readlines()
    f.close()
    n_lines = len(lines)               # number of lines in a file

    # read header line
    line = lines[0]
    param_names = line.strip().split('|')[0].split()
    qoi_names   = line.strip().split('|')[1].split()
    all_names   = param_names + qoi_names
    
    # read simulations
    self.all_sims = []                      # simulation
    for i_line in range(1,n_lines):
        
      # parse parameter values
      params = lines[i_line].strip().split('|')[0].split()
      for i, param in enumerate(params):
        if i == 0:
          # i == 0 is the simulation number
          params[i] = int(param)
        else:
          params[i] = float(param)
          
      # parse qoi values
      qois = lines[i_line].strip().split('|')[1].split()
      for i, qoi in enumerate(qois):
        qois[i] = float(qoi)
        
      self.all_sims.append(params+qois)
        
    self.filename_in = fname_in
    self.param_names = param_names
    self.qoi_names = qoi_names
    self.all_names = all_names
    self.np_all_sims = np.array(self.all_sims)
    
    self.__get_names_for_error_types()

  def set_qoi_type(self,qoi_type):
    self.qoi_type = qoi_type
    self.qois = []
    if qoi_type == 'abs_err':
        self.qois = self.names_abserr
    elif qoi_type == "nabs_err":
        self.qois = self.names_nabserr
    elif qoi_type == "sqerr":
        self.qois = self.names_sqerr
    elif qoi_type == "nsqerr":
        self.qois = self.names_nsqerr
    else:
        raise ValueError("unknown_qoi_type")
  
  def create_dataset_for_pareto_analysis(self):
    self.dataset = []
    for n, sim in enumerate(self.all_sims):
        self.dataset.append(pyflamestk.pareto.Datapoint(n))
        for qoi in self.qois:
          idx = self.all_names.index(qoi)
          self.dataset[n].addNumber(-sim[idx])
          
  def calculate_pareto_set(self):
    pyflamestk.pareto.bruteforce_algo(self.dataset)
    self.pareto_set_ids = []
    self.pareto_set = []
    for s in self.dataset:
      if s.paretoStatus == 1:
        self.pareto_set_ids.append(s.id)
        self.pareto_set.append(s.vec)
    self.pareto_set = -np.array(self.pareto_set)
    self.np_pareto_set = self.np_all_sims[self.pareto_set_ids]
    self.pareto_mean = np.mean(self.np_pareto_set)
    self.pareto_cov  = np.cov(self.np_pareto_set)
    write_pareto_set(param_names=self.param_names,
                     qoi_names=self.param_names,
                     pareto_set=self.pareto_set,
                     fname_out=self.fname_pareto_out)

  def add_performance_constraint(self,metric_name,metric_value):
    self.performance_requirements[metric_name] = metric_value

  def apply_performance_constraints(self):
    # start with the full pareto set and then remove elements which do 
    # not meat the performane criteria
    n_sims, n_qoi = self.np_pareto_set.shape
    self.np_pareto_set_cull = np.copy(self.np_pareto_set)

    #determine which rows to delete
    rows_to_delete = []
    for idx in range(n_sims):
      ps = self.np_pareto_set_cull[idx,:]
      is_delete_row = False
      for qoi_name in self.performance_requirements.keys():
        qoi_idx = self.all_names.index(qoi_name)
        if ps[qoi_idx] > self.performance_requirements[qoi_name]:
          is_delete_row = True 
      if is_delete_row:
        rows_to_delete.append(idx)
        
    # remove rows which do not meet performance criteria
    self.np_pareto_set_cull = np.delete(self.np_pareto_set_cull,
                                        rows_to_delete,
                                        axis=0)
    self.pareto_cull_mean = np.mean(self.np_pareto_set_cull)
    self.pareto_cull_cov  = np.cov(self.np_pareto_set_cull)

  def create_all_histograms(self):
    for name in self.param_names:
      self.create_histogram(name)

  def create_histogram(self,name):
    plt.figure()
    x_idx = self.all_names.index(name)
    plt.hist(self.np_all_sims[:,x_idx])
    plt.show()

  def create_all_pareto_plots(self,qoi_list):
    for i, qoi_name_i in enumerate(qoi_list):         
      for j, qoi_name_j in enumerate(qoi_list):
        if i < j and i != j:
          print("{} {}".format(qoi_name_i,qoi_name_j))
          x_label = qoi_name_i
          y_label = qoi_name_j
          x_idx   = self.all_names.index(x_label)
          y_idx   = self.all_names.index(y_label)
          pareto_front = pareto_frontier_2d(self.np_all_sims[:,x_idx],
                                            self.np_all_sims[:,y_idx],
                                            maxX = False, maxY= False)
          fig, ax = plt.subplots()
          ax.scatter(self.np_all_sims[:,x_idx],
                     self.np_all_sims[:,y_idx],
                     label='dominated')
          ax.scatter(self.pareto_set[:,self.qois.index(x_label)],
                     self.pareto_set[:,self.qois.index(y_label)],label = 'pareto', color='y',)
          ax.plot(pareto_front[0],
                  pareto_front[1],
                  color='r',
                  linewidth=2)
          legend = ax.legend(loc="upper right")
          plt.axis([min(pareto_front[0]),
                    max(pareto_front[0]),
                    min(pareto_front[1]),
                    max(pareto_front[1])])
          plt.xlabel(x_label)
          plt.ylabel(y_label)
          plt.show()

      
  def __get_names_for_error_types(self):    
    # get names for different types
    self.names_abserr = []  # initialize, array for q_est - q
    self.names_nabserr = [] # initialize, array for (q_est - q)/q
    self.names_sqerr  = []  # initialize, array for q_est - q
    self.names_nsqerr = []  # initialize, array for (q_est -q)/q
    
    for name in self.all_names:
      if name.endswith("_abserr"):
        self.names_abserr.append(name)
      if name.endswith("_nabserr"):
        self.names_nabserr.append(name)
      if name.endswith("_sqerr"):
        self.names_sqerr.append(name)
      if name.endswith("_nsqerr"):
        self.names_nsqerr.append(name)

#------------------------------------------------------------------------------
# %%SCRIPT STARTS HERE

fname_results_1 = "results_10000b.out"
fname_results_2 = "results_10000b.out"
qoi_type        = "abs_err"
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

sim_results = SimulationResults()
sim_results.read_simulation_results(fname_in=fname_results_1)
sim_results.set_qoi_type(qoi_type)
sim_results.create_dataset_for_pareto_analysis()
sim_results.calculate_pareto_set()

sim_results.add_performance_constraint('MgO_NaCl_latt_abserr', 0.05)
sim_results.add_performance_constraint('MgO_NaCl_B_abserr', 50.)
sim_results.add_performance_constraint('MgO_NaCl_G_abserr', 50.)
sim_results.apply_performance_constraints()

#%%
## read in lines
#f = open(fname_results_1)
#lines = f.readlines()
#f.close()
#
#n_lines  = len(lines)     # number of lines in a file
#
## read header line
#line = lines[0]
#param_names = line.strip().split('|')[0].split()
#qoi_names   = line.strip().split('|')[1].split()
#all_names   = param_names + qoi_names
#
## read simulations
#all_sims = []   #initialize
#for i_line in range(1,n_lines):
#    
#  params = lines[i_line].strip().split('|')[0].split()
#  for i, param in enumerate(params):
#      if i == 0:
#          params[i] = int(param)
#      else:
#          params[i] = float(param)
#          
#  qois   = lines[i_line].strip().split('|')[1].split()
#  for i, qoi in enumerate(qois):
#      qois[i] = float(qoi)
#  all_sims.append(params+qois)
#np_all_sims = np.array(all_sims)
#
## get names for different types
#names_abserr = []  # initialize, array for q_est - q
#names_nabserr = [] # initialize, array for (q_est - q)/q
#names_sqerr  = []  # initialize, array for q_est - q
#names_nsqerr = []  # initialize, array for (q_est -q)/q
#for name in all_names:
#    if name.endswith("_abserr"):
#        names_abserr.append(name)
#    if name.endswith("_nabserr"):
#        names_nabserr.append(name)
#    if name.endswith("_sqerr"):
#        names_sqerr.append(name)
#    if name.endswith("_nsqerr"):
#        names_nsqerr.append(name)
#
#qois = []
#if qoi_type == 'abs_err':
#    qois = names_abserr
#elif qoi_type == "nabs_err":
#    qois = names_nabserr
#elif qoi_type == "sqerr":
#    qois = names_sqerr
#elif qoi_type == "nsqerr":
#    qois = names_nsqerr
#else:
#    raise ValueError("unknown_qoi_type")

#------------------------------------------------------------------------------
# create dataset for pareto analysis    
#dataset = []
#for n, sim in enumerate(all_sims):
#    dataset.append(pyflamestk.pareto.Datapoint(n))
#    for qoi in qois:
#      idx = all_names.index(qoi)
#      dataset[n].addNumber(-sim[idx])
#
##------------------------------------------------------------------------------
## calculate pareto set
#pyflamestk.pareto.bruteforce_algo(dataset)
#pareto_set_ids = []
#pareto_set = []
#for s in dataset:
#    if s.paretoStatus == 1:
#        pareto_set_ids.append(s.id)
#        pareto_set.append(s.vec)
#pareto_set = -np.array(pareto_set)
#write_pareto_set(param_names=param_names,
#                 qoi_names=param_names,
#                 pareto_set=pareto_set,
#                 fname_out='pareto.out')
                 
#--------------------------------------------------------------------
# filter pareto set by performance requirements
    

#%%
#------------------------------------------------------------------------------
# make plots

# parameter histograms
#for name in param_names:
#    create_histogram(name)
    
# simulations are arranged in array[row,columns]
    
#%% CREATE PARETO CURVE GRAPHS
#qoi_list = ['MgO_NaCl_latt_abserr',
#            'MgO_NaCl_c11_abserr',
#            'MgO_NaCl_c12_abserr',
#            'MgO_NaCl_c44_abserr',
#            'MgO_NaCl_B_abserr',
#            'MgO_NaCl_G_abserr']
#
#for i, qoi_name_i in enumerate(qoi_list):         
#  for j, qoi_name_j in enumerate(qoi_list):
#    if i < j and i != j:
#      print("{} {}".format(qoi_name_i,qoi_name_j))
#      x_label = qoi_name_i
#      y_label = qoi_name_j
#      x_idx   = all_names.index(x_label)
#      y_idx   = all_names.index(y_label)
#      pareto_front = pareto_frontier_2d(np_all_sims[:,x_idx],
#                                        np_all_sims[:,y_idx],
#                                        maxX = False, maxY= False)
#      fig, ax = plt.subplots()
#      ax.scatter(np_all_sims[:,x_idx],
#                 np_all_sims[:,y_idx],
#                 label='dominated')
#      ax.scatter(pareto_set[:,qois.index(x_label)],
#                 pareto_set[:,qois.index(y_label)],label = 'pareto', color='y',)
#      ax.plot(pareto_front[0],
#              pareto_front[1],
#              color='r',
#              linewidth=2)
#      legend = ax.legend(loc="upper right")
#      plt.axis([min(pareto_front[0]),
#                max(pareto_front[0]),
#                min(pareto_front[1]),
#                max(pareto_front[1])])
#      plt.xlabel(x_label)
#      plt.ylabel(y_label)
#      plt.show()
#%% MAKE 
n = len(qois)
for i in range(n):
    for j in range(n):
        if i < j and i != j:
            x_label = qois[i]
            y_label = qois[j]
            x_idx   = all_names.index(x_label)
            y_idx   = all_names.index(y_label)
            pareto_front = pareto_frontier_2d(np_all_sims[:,x_idx],
                                              np_all_sims[:,y_idx],
                                              maxX = False, maxY= False)
            plt.figure()
            plt.scatter(np_all_sims[:,x_idx],
                        np_all_sims[:,y_idx])
            plt.scatter(pareto_set[:,qois.index(x_label)],
                        pareto_set[:,qois.index(y_label)],color='y')
            plt.plot(pareto_front[0],pareto_front[1])
            plt.axis([min(pareto_front[0]),
                      max(pareto_front[0]),
                      min(pareto_front[1]),
                      max(pareto_front[1])])
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.show()

#%% HISTOGRAMS
param_names_id = []
for i in param_names:
    param_names_id.append(all_names.index(i))
param_names_id.remove(param_names_id[0])
param_names_id.remove(param_names_id[3])

import scipy.stats

np_pareto_set = np_all_sims[pareto_set_ids]
for name in param_names:
    plt.figure()
    plt.subplot(211)
    x_idx = all_names.index(name)
    plt.hist(np_pareto_set[:,x_idx],bins=40)

    plt.subplot(212)
    try:
      data = np_pareto_set[:,x_idx]
      density = scipy.stats.gaussian_kde(data)
      x = np.linspace(min(data), max(data), 1000)
      
      #density.covariance_factor = lambda : .25
      #density._compute_covariance()
      plt.plot(x, density(x))
      plt.xlabel(name)

    except np.linalg.linalg.LinAlgError:
       print("error") 
    plt.show()

#%%
    
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
  #ax.plot(m1, m2, 'k.', markersize=2)
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  #plt.set_xlim([xmin, xmax])
  #plt.set_ylim([ymin, ymax])
  cb = plt.colorbar()
  cb.set_label("freq")
  plt.show()

param_list = ['chrg_Mg',  'chrg_O',     'p_MgMg_a',
              'p_MgMg_c', 'p_MgMg_rho', 'p_OO_a',
              'p_OO_c',   'p_OO_rho',   'p_MgO_a',
              'p_MgO_c',  'p_MgO_rho']
              
for pn_i, x_label in enumerate(param_list):
    for pn_j, y_label in enumerate(param_list):
        if pn_i < pn_j and pn_i != pn_j:
            x_data  = np_pareto_set[:,all_names.index(x_label)]
            y_data  = np_pareto_set[:,all_names.index(y_label)]
            make_density_plot(x_data,y_data,x_label,y_label)
            try:
              make_jpdf_plot(x_data,y_data,x_label,y_label)
            except np.linalg.linalg.LinAlgError:
              print("Kernel Density Estimator produced LinAlgErr")


#%% Cull the pareto set by the performance requirements
performance_requirements = {}
performance_requirements['MgO_NaCl_latt_abserr'] = 0.05
performance_requirements['MgO_NaCl_B_abserr'] = 50
performance_requirements['MgO_NaCl_G_abserr'] = 50

n_sims, n_qoi = np_pareto_set.shape
np_pareto_set_cull = np.copy(np_pareto_set)
rows_to_delete = []
for idx in range(n_sims):
  ps = np_pareto_set_cull[idx,:]
  is_delete_row = False
  for qoi_name in performance_requirements.keys():
    qoi_idx = all_names.index(qoi_name)
    if ps[qoi_idx] > performance_requirements[qoi_name]:
      is_delete_row = True 
  if is_delete_row:
    rows_to_delete.append(idx)
    
np_pareto_set_cull = np.delete(np_pareto_set_cull,rows_to_delete,axis=0)

param_list = ['chrg_Mg',  'chrg_O',     'p_MgMg_a',
              'p_MgMg_c', 'p_MgMg_rho', 'p_OO_a',
              'p_OO_c',   'p_OO_rho',   'p_MgO_a',
              'p_MgO_c',  'p_MgO_rho']
              
for pn_i, x_label in enumerate(param_list):
    for pn_j, y_label in enumerate(param_list):
        if pn_i < pn_j and pn_i != pn_j:
            x_data  = np_pareto_set_cull[:,all_names.index(x_label)]
            y_data  = np_pareto_set_cull[:,all_names.index(y_label)]
            make_density_plot(x_data,y_data,x_label,y_label)
            try:
              make_jpdf_plot(x_data,y_data,x_label,y_label)
            except np.linalg.linalg.LinAlgError:
              print("Kernel Density Estimator produced LinAlgErr")
              
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
plt.violinplot(pareto_set[:,0], showmeans=False, showmedians=True)
plt.violinplot(pareto_set[:,1], showmeans=False, showmedians=True)
plt.violinplot(pareto_set[:,2], showmeans=False, showmedians=True)
plt.violinplot(pareto_set[:,3], showmeans=False, showmedians=True)
plt.violinplot(pareto_set[:,4], showmeans=False, showmedians=True)
plt.violinplot(pareto_set[:,5], showmeans=False, showmedians=True)

min_error = np.zeros((1,6))
min_error[:,1] = 1
min_error[:,2] = 1
min_error[:,3] = 1