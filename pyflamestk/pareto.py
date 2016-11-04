#!/usr/bin/env python
import pyflamestk.pyposmat
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats
import numpy as np
import sys, getopt
import copy

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 

        original code: http://stackoverflow.com/questions/11882393/matplotlib-disregard-outliers-when-plotting
        date: 10/19/2016
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
    
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

    # write header line    
    str_out = ""
    for name in param_names:
        str_out += "{} ".format(name)
    str_out += "| "
    for name in qoi_names:
        str_out += "{} ".format(name)
    str_out += "\n"
    f.write(str_out)

    # write body    
    for myset in pareto_set:
        str_out = ""
        for item in myset:
            str_out += "{} ".format(item)
        str_out += "\n"
    f.write(str_out)
    
    f.close()
    
class ParameterFileReader:
    """
    Args:
        fname_param_in (str): filename of the parameter file
        
    Attributes:
        param_names     (python array) parameter names.  same index as params.
        params          (numpy.array[sim_id,param_id])
    """
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
            
class SimulationResults:
  """
  
  Args:
      n_simulations (int): number of simulations read from the output file
      qoi_type (str): supported qoi types are (1) 'abserr' - absolute error,
                      (2) 'nabserr' - normalized absolute error, 
                      (3)'sqerr' - square error, and (4) 'nsqerr' - normalized 
                      square error
  """
  def __init__(self):      
    self.performance_requirements = {}
    self.fname_pareto_out = "pareto.out"

    # TODO: Initialization of variables here is very sloppy
    #       encapsulate variables that can be encapsulated

    # initialize variables
    self.all_sims    = []     # python array of all simulation data
    self.pareto_set  = []     # python array of the pareto set


    # initialize variables [ATTRIBUTES]
    self.n_simulations = 0
    self.qoi_type = 'abserr'  # supported types: abserr, nabserr, sqerr, nsqerr
    self.qoi_names   = []     # array of qoi names
    self.param_names = []     # array of parameter names
    self.qois = []
    self.qoi_keys = []        # TODO: horrible name fix
    self.np_pareto_set = []
    self.np_pareto_set_ids = []
    self.np_all_sims = []     # numpy array of all simulation data
    self.n_resamples = 0

    self.performance_requirements = {}
    
    self.pareto_dataset = []
    self.pareto_set = []
    self.pareto_set_id = []

  def read_simulation_results(self,fname_in):
      """read simulations results from a file into a memory.
      
      Args:
          fname_in (str): the filename containing the simulation results from
                          LAMMPS simulations
      """

      self.filename_in = fname_in
      
      # read file
      f = open(self.filename_in)
      lines = f.readlines()
      f.close()
      n_lines = len(lines)               # number of lines in a file
      self.n_lines = len(lines)
      self.n_simulations = n_lines - 1
      self.n_resamples = n_lines - 1
    
      # read header line
      line = lines[0]
      param_names = line.strip().split('|')[0].split()
      qoi_names   = line.strip().split('|')[1].split()
      all_names   = param_names + qoi_names
    
      # read simulations
      self.all_sims = [] # initialization                     
      for i_line in range(1,n_lines):
        
          # the '|' separates the parameter values from the qois
          params = lines[i_line].strip().split('|')[0].split()
          qois = lines[i_line].strip().split('|')[1].split()
          
          # parse parameters
          for i, param in enumerate(params):
              if i == 0:
                  params[i] = int(param)
              else:
                  params[i] = float(param)
            
          # parse qoi values
          for i, qoi in enumerate(qois):
              qois[i] = float(qoi)
        
          self.all_sims.append(params+qois)
        
      self.param_names = param_names
      self.qoi_names   = qoi_names
      self.all_names   = all_names
      self.np_all_sims = np.array(self.all_sims)
    
      self.__get_names_for_error_types()

  def set_qoi_type(self,qoi_type):
    self.qoi_type = qoi_type
    self.qois = []
    if qoi_type == 'abserr':
        self.qois = self.names_abserr
    elif qoi_type == "nabserr":
        self.qois = self.names_nabserr
    elif qoi_type == "sqerr":
        self.qois = self.names_sqerr
    elif qoi_type == "nsqerr":
        self.qois = self.names_nsqerr
    else:
        raise ValueError("unknown_qoi_type")
  
  def __create_dataset_for_pareto_analysis(self, qoi_keys=""):
      """
      qoi_keys      If not set, will use all qois
      """

      # Pareto set for all qois
      if qoi_keys == "":
          self.qoi_keys = copy.deepcopy(self.qois)
      else:
          self.qoi_keys = copy.deepcopy(qoi_keys)
          for i in range(len(self.qoi_keys)):
              if self.qoi_type == 'abserr':
                  self.qoi_keys[i] = "{}_abserr".format(self.qoi_keys[i])
              elif self.qoi_type == 'nabserr':
                  self.qoi_keys[i] = "{}_nabserr".format(self.qoi_keys[i])
              elif self.qoi_type == 'sqerr':
                  self.qoi_keys[i] = "{}_sqerr".format(self.qoi_keys[i])
              elif self.qoi_type == 'nsqerr':
                  self.qoi_keys[i] = "{}_nsqerr".format(self.qoi_keys[i])
              else:
                  raise ValueError("unknown_qoi_type")
                  
      # check to see if each key exists
      for qoi_key in self.qoi_keys:
          if not qoi_key in self.qois:
              errmsg = "qoi_key, {}, does not exist"
              errmsg = errmsg.format(qoi_key)
              raise ValueError(errmsg)
              
      # make dataset 
      self.pareto_dataset = []
      for n, sim in enumerate(self.all_sims):
          self.pareto_dataset.append(pyflamestk.pareto.Datapoint(n))
          for qoi in self.qoi_keys:
              idx = self.all_names.index(qoi)
              self.pareto_dataset[n].addNumber(-sim[idx])
         
  def calculate_pareto_set(self,qoi_keys):
      my_qoi_keys = copy.deepcopy(qoi_keys)

      self.__create_dataset_for_pareto_analysis(qoi_keys=my_qoi_keys)
      print("calculating pareto set...")
      pyflamestk.pareto.bruteforce_algo(self.pareto_dataset)
      self.pareto_set_ids = []
      self.pareto_set = []
    
      # mark pareto set
      for s in self.pareto_dataset:
        if s.paretoStatus == 1:
          self.pareto_set_ids.append(s.id)
          self.pareto_set.append(s.vec)
      self.pareto_set = -np.array(self.pareto_set)
      self.np_pareto_set = self.np_all_sims[self.pareto_set_ids]
      self.pareto_mean = np.mean(self.np_pareto_set)
      self.pareto_cov  = np.cov(self.np_pareto_set)

  def calculate_parameter_estimates(self,param_list):
      params = copy.deepcopy(param_list)
      self.param_estimates = {}
      for param in params:
          self.param_estimates[param] = {}
          self.param_estimates[param]['all'] = {}
          self.param_estimates[param]['pareto'] = {}
          self.param_estimates[param]['pareto_cull'] = {}
  
  def calculate_qoi_estimates(self,qoi_keys):
      qois = copy.deepcopy(qoi_keys)
      set_types = ['all','pareto','pareto_cull']
      self.qoi_estimates = {}
      for qoi in qois:
          self.qoi_estimates[qoi] = {}
          for set_type in set_types:
              #TODO
              mean = 0
              std  = 1
              self.qoi_estimates[qoi][set_type] = {}
              self.qoi_estimates[qoi][set_type]['mean'] = mean
              self.qoi_estimates[qoi][set_type]['std'] = std
          
  def cull_by_percentile(self,pct_kept=10.):
        """
        
        Arguments:
        pct_kept (float, 10.0) - number between 1 and 100 indicating the pct of simulations
            within the Pareto set which should be kept
                        
        Returns:
        
        a numpy array with observations indexed in rows, and parameters
        and quantities of interst indexed in columns.  The column index is the
        same as the array in "self.all_names".
        
        Note:
        
        TODO:
        
        1)   A Newton-Ralphson method to get more accurate performance requirements
        to prevent over culling of the Pareto set.
        """
        
        if 0 <= pct_kept <= 100.:
            errmsg = "pct_kept must be between 1 and 100, the value {} was passed."
            errmsg = errmsg.format(pct_kept)
            raise ValueError(errmsg)
        else:
            self.pct_kept = pct_kept
            
            
        qoi_keys = self.qoi_keys
        self.performance_requirements = {}
        for qoi_key in qoi_keys:
            self.performance_requirements[qoi_key] = 0.
        n_sims, n_qoi = self.np_pareto_set.shape        
    
        # intialize variables
        pctl_threshold = 100        # searching for 100% within the Pareto set
                                    # to 0% in the pareto set
        is_culled = False           # intialize
        
        while not is_culled:
            rows_to_delete = []
            pctl_threshold -= 1
            # calculate percentile cutoffs
            for qoi_key in self.performance_requirements.keys():
                if pctl_threshold < 0:
                    errmsg = "While searching for the pctl_threshold, the percentile error dropped below zero resulting in an error."
                    raise ValueError(errmsg)
                    # TODO: remove
                    # print(msg_out)
                    # print("percentile must be greater than 0")
                    # is_culled = True
                else:
                    qoi_data = self.get_data_by_name(qoi_key, 'pareto')
                    #TODO: remove
                    #qoi_id = self.all_names.index(qoi_key)
                    # qoi_data = self.np_pareto_set[:,qoi_id]
                    self.performance_requirements[qoi_key] = np.percentile(qoi_data, pctl_threshold)

            # cull the pareto set by the performance requirements
            for qoi_key in self.performance_requirements.keys():        
                np_pareto_set_cull = np.copy(self.np_pareto_set)
                for idx in range(n_sims):
                    ps = np_pareto_set_cull[idx,:]
                    is_delete_row = False
                    for qoi_name in self.performance_requirements.keys():
                        qoi_idx = self.all_names.index(qoi_name)
                        if ps[qoi_idx] > self.performance_requirements[qoi_name]:
                            is_delete_row = True 
                    if is_delete_row:
                        rows_to_delete.append(idx)
            
            # check to see if the pareto set has been sufficiently culled
            n_culled = len(rows_to_delete)
            pct_culled = float(n_culled/n_sims)
            if pct_kept/100. > 1 - pct_culled:
                is_culled = True
        
        # delete rows not in the pareto set
        self.np_pareto_set_cull = np.copy(self.np_pareto_set)
        self.np_pareto_set_cull = np.delete(np_pareto_set_cull,rows_to_delete,axis=0)

        # TODO: this should be moved to a reporting section
        msg_out = ""                # intialize output string
        msg_out += "n_pareto_set        = {}\n".format(n_sims)
        msg_out += "n_rows_culled       = {}\n".format(n_culled)
        msg_out += "n_culled_pareto_set = {}\n".format(n_sims - n_culled)
        msg_out += "performance constraints:\n"
        for qoi_key in self.performance_requirements.keys():
            msg_out += "\t {} < {:0.4f}\n".format(qoi_key,self.performance_requirements[qoi_key])    
        print(msg_out)
        
        return self.np_pareto_set_cull.copy()

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
    
  def get_data_by_name(self,name,ds_type):
      """
      Arguments:
      
      name (str) - string of parameter or quantity of interest
      ds_type (str) - string of which dataset we are taking the data from.
          The ds_types which are supported are: all, pareto, pareto_culled
          
      Returns:
      
      a numpy array of the data asked for
          
      """
      
      idx = self.all_names.index(name)
      
      if ds_type == 'all':
          return self.np_all_sims[:,idx]
      elif ds_type == 'pareto':
          return self.np_pareto_set[:,idx]
      elif ds_type == 'pareto_cull':
          return self.np_pareto_set_cull[:,idx]
      else:
          errmsg = "ds_type must be either all, pareto, or pareto_cull"
          raise ValueError(errmsg)
 
  def create_2d_pareto_plot(self,
                            qoi_name_1,
                            qoi_name_2,
                            error_type = "",
                            show_dominated = True,
                            show_pareto = True, 
                            show_culled = True,
                            show_2d_pareto_curve = True,
                            show_2d_culled_curve = True):   
        """
        Arguments:
        
        qoi_name_1 (string) - the name of the first quantity of interest to be
            plotted on the x-axis
        qoi_name_2 (string) - the name of the second quantity of interest to be
            plotted on the y-axis
        error_type (string: "" ) - not implemented yet
        show_dominated (bool, True) - shows dominated points when set to true
        show_pareto (bool,True) - shows Pareto points when set to true.
        show_culled (bool,True) - shows the remaining Pareto points when the 
            Pareto set is culled by by performance requirements.
        show_2d_pareto_curve (bool,True) - shows a 2d slice of the Pareto curve
            when set to True.
        show_2d_pareto_cuve (bool,True) - shows a 2d slice of the culled Pareto
            curve when set to false.
            
        Returns:
        
        Nothing
            
        Notes:
        
        When the dominated points are plotted, it actually plots all the points
        in the dataset.  However, these points are then plotted over by the
        Pareto points.
        """
        
        
        x_label = qoi_name_1
        y_label = qoi_name_2
                                   
        x_data_all = self.get_data_by_name(x_label, 'all')
        y_data_all = self.get_data_by_name(y_label, 'all')
        
        if show_pareto == True:
            x_data_pareto = self.get_data_by_name(x_label, 'pareto')
            y_data_pareto = self.get_data_by_name(y_label, 'pareto')

            if show_2d_pareto_curve == True:
                pareto_2d = pyflamestk.pareto.pareto_frontier_2d(x_data_all,
                                                                 y_data_all,
                                                                 maxX = False, maxY= False) 
            
        if show_culled == True:
            x_data_cull = self.get_data_by_name(x_label, 'pareto_cull')
            y_data_cull = self.get_data_by_name(y_label, 'pareto_cull')

            if show_2d_pareto_curve == True:
                cull_2d = pyflamestk.pareto.pareto_frontier_2d(x_data_cull,
                                                               y_data_cull,
                                                               maxX = False, maxY= False)         
        # plot results
        plt.figure()
    
        if show_dominated == True:
            c = 'b'
            plt.scatter(x_data_all, y_data_all, color = c)
            
        if show_pareto == True:
            c = 'y'
            plt.scatter(x_data_pareto,y_data_pareto, color=c)
            if show_2d_pareto_curve == True:
                plt.plot(pareto_2d[0],pareto_2d[1],color=c)

        if show_culled == True:
            c = 'g'
            plt.scatter(x_data_cull,  y_data_cull, color = c)
            if show_2d_pareto_curve == True:
                plt.plot(cull_2d[0],cull_2d[1],color=c)

        # determine axis
        xmin = min(pareto_2d[0])
        xmax = max(pareto_2d[0])
        ymin = min(pareto_2d[1])
        ymax = max(pareto_2d[1])
        plt.axis([xmin,xmax,ymin,ymax])
        
        # add labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        
        # display graph
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
    
class ParetoSet:
  def __init__(self):
    pass  

class Datapoint:
    """Defines a point in K-dimensional space"""
    def __init__(self,id):
        self.id = id # datapoint id (0,..N-1)
        self.vec = [] # the K-dim vector
        self.paretoStatus = -1 # -1=dont know, 1=pareto, 0=not pareto
        self.dominatedCount = 0 # number of datapoints that dominate this point
        self.dominatingSet = [] # set of vectors this one is dominating

    def addNumber(self,num):
        """Adds a number to one dimension of this datapoint"""
        self.vec.append(num)

    def addToDominatingSet(self,id2):
        """Add id of of dominating point"""
        self.dominatingSet.append(id2)

    def dominates(self,other):
        """Returns true if self[k]>=other[k] for all k and self[k]>other[k] for at least one k"""
        assert isinstance(other,Datapoint)
        gte=0 # count of self[k]>=other[k]
        gt=0 # count of self[k]>other[k]
        for k in range(len(self.vec)):
            if self.vec[k] >= other.vec[k]:
                gte+=1
                if self.vec[k] > other.vec[k]:
                    gt+=1
            
        return (gte==len(self.vec) and (gt>0))

    def __repr__(self):
        return self.vec.__repr__()+": "+str(self.paretoStatus)

def bruteforce_algo(dataset):
    num_pareto = 0
    
    # pairwise comparisons
    for n in range(len(dataset)):
        if np.mod(n,100) == 0:
          print("\t n={}".format(n))
        for m in range(len(dataset)):
            if dataset[m].dominates(dataset[n]):
                dataset[n].dominatedCount+=1
                dataset[m].addToDominatingSet(n)

    # find first pareto front
    for n in range(len(dataset)):
        if dataset[n].dominatedCount == 0:
            dataset[n].paretoStatus = 1
            num_pareto += 1
        else:
            dataset[n].paretoStatus = 0
                
def nondominated_sort(dataset):
    """Nondominated Sorting, generates ranking w/ higher number = better pareto front"""
    numPareto = 0

    # pairwise comparisons
    for n in range(len(dataset)):
        for m in range(len(dataset)):
            if dataset[m].dominates(dataset[n]):
                dataset[n].dominatedCount+=1
                dataset[m].addToDominatingSet(n)

    # find first pareto front
    front = []
    front2 = []
    tmpLevel = -10 # temporary value for Pareto level, will re-adjust later
    for n in range(len(dataset)):
        if dataset[n].dominatedCount == 0:
            dataset[n].paretoStatus = tmpLevel
            front.append(n)
            numPareto+=1

    # iteratively peel off pareto fronts
    while len(front) != 0:
        tmpLevel-=1
        for f in front:
            for s in dataset[f].dominatingSet:
                dataset[s].dominatedCount -= 1
                if dataset[s].dominatedCount == 0:
                    front2.append(s)
                    dataset[s].paretoStatus = tmpLevel
        front = front2 
        front2 = []

    # re-adjust pareto level
    for n in range(len(dataset)):
        oldLevel = dataset[n].paretoStatus
        if oldLevel != -1:
            dataset[n].paretoStatus = oldLevel-tmpLevel-1

    return numPareto


def create_dataset(raw_vectors):
    """Given a list of vectors, create list of datapoints"""
    dataset = []
    for k in range(len(raw_vectors)):
        for n,v in enumerate(raw_vectors[k]):
            if k == 0:
                dataset.append(Datapoint(n))
            dataset[n].addNumber(v)
    return dataset


def readfile(filename,multiplier=1.0):
    """Reads a vector file (objective values in one dimension)"""
    with open(filename,'r') as f:
        lines = f.readlines()
    vec = [multiplier*float(a.strip()) for a in lines]
    return vec

def make_histograms_params_combined(sim_results,params):
    #TODO: Implement this function
    raise NotImplementedError("this function has not been implemented yet")
    #print("making combined histograms for parameters")
    #n = len(params)
    #for i in range(n):
    #    param_label = params[i]
        
        
def make_histograms_qois_combined(sim_results,qois):
    print("making combined histograms")
    n = len(qois)
    for i in range(n):
        qoi_label = qois[i]
        qoi_value = sim_results.qoi_values[qoi_label]
        qoi_idx = sim_results.all_names.index(qoi_label)

        # all simulated points
        qoi_values_all = np.copy(sim_results.np_all_sims[:,qoi_idx])
        qoi_values_all = qoi_values_all[~pyflamestk.pareto.is_outlier(qoi_values_all)]
        qoi_mean_all = np.mean(qoi_values_all)
        qoi_std_all  = np.std(qoi_values_all)

        # points in the pareto set
        qoi_values_pareto = np.copy(sim_results.np_pareto_set[:,qoi_idx])
        qoi_values_pareto = qoi_values_pareto[~pyflamestk.pareto.is_outlier(qoi_values_pareto)]
        qoi_mean_pareto = np.mean(qoi_values_pareto)
        qoi_std_pareto = np.std(qoi_values_pareto)

        print("{}, mean: {}, std:{}".format(qoi_label,qoi_mean_all,qoi_std_all))
        fname = "hist_qoi_{}.jpg".format(qoi_label)
        print("filename: {}".format(fname))
        plt.figure()
        # this plot has two superimposed histograms, the first histogram
        # contains all simulated predictions, the second histogram only has
        # predictions which are part of the pareto set
        fig, ax1 = plt.subplots()
        common_params = dict(bins=20,
                             range=(min(qoi_values_all),max(qoi_values_all)),
                             alpha=1.0)
        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_all,    facecolor='g', **common_params)
        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_pareto, facecolor='y', **common_params)
        plt.xlabel(qoi_label)
        plt.axvline(x=qoi_value,color='black',linewidth=2)
        plt.ylabel('frequency')

        # add density function
        ax2 = ax1.twinx()   
        ax2.plot(qoi_bins,
                 mlab.normpdf(qoi_bins,qoi_mean_all,qoi_std_all),
                 'r--', linewidth=2)
        ax2.set_ylabel('probability density')
        plt.show()

        fname = "hist_qoi_{}.jpg".format(qoi_label)
        print("{}, mean: {}, std:{}".format(qoi_label,qoi_mean_pareto,qoi_std_pareto))
        print("filename = {}".format(fname))
        plt.figure()
        # this plot only contains one histogram contain the predictions which
        # are part of the pareto set
        fig, ax1 = plt.subplots()

        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_pareto, bins=20, facecolor='y', alpha=1.0)
        plt.axvline(x=qoi_value,color='black',linewidth=2)
        plt.xlabel(qoi_label)
        plt.ylabel('frequency')
        
        # add density function
        ax2 = ax1.twinx()   
        ax2.plot(qoi_bins,
                 mlab.normpdf(qoi_bins,qoi_mean_pareto,qoi_std_pareto),
                 'r--', linewidth=2)
        ax2.set_ylabel('probability density')

        plt.show()
        
def make_histograms_parameters_combined(sim_results):
    print("making parameter distribution histograms")
    
    # create dictionary item for parameter names and parameter ids
    param_names_id = {}
    for idx, param_name in enumerate(sim_results.param_names):
        if param_name != 'sim_id':
            param_names_id[param_name] = sim_results.all_names.index(param_name)
            
    # iterate over the dictionary
    for param_name, param_id in param_names_id.items():

        param_values_all = sim_results.np_all_sims[:,param_id]
        param_all_mu  = param_values_all.mean()
        param_all_std = param_values_all.std()

        param_values_pareto = sim_results.np_pareto_set[:,param_id]
        param_pareto_mu  = param_values_pareto.mean()
        param_pareto_std = param_values_pareto.std()
        
        param_values_cull = sim_results.np_pareto_set_cull[:,param_id]
        param_cull_mu  = param_values_cull.mean()
        param_cull_std = param_values_cull.std()        

        print("all:    {} {} {} {}".format(param_name, param_id, param_all_mu, param_all_std))
        print("pareto: {} {} {} {}".format(param_name, param_id, param_pareto_mu, param_pareto_std))
        print("cull:   {} {} {} {}".format(param_name, param_id, param_cull_mu, param_cull_std))

        plt.figure()
        plt.hist(param_values_all, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('all simulations')
        plt.show()
        
        plt.figure()
        plt.hist(param_values_pareto, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('pareto set')
        plt.show()

        plt.figure()
        plt.hist(param_values_cull, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('pareto set with performance constraints')
        plt.show()
        
def make_jpdf_plot(x_data,
                   y_data,
                   x_label,
                   y_label, 
                   axis="", 
                   title=""):
    """ Make a joint probability density plot (jpdf) plot using a kernel
    density estimate.
    
    Arguments:
    
    x_data (numpy.array)
    y_data (numpy.array)
    x_label (string)
    y_label (string)
    axis (array)
    title (string)
    """
    
    xmin = 0.
    ymax = 0.
    ymin = 0.
    ymax = 0.
    if axis == "":
        xmin = x_data.min()
        xmax = x_data.max()
        ymin = y_data.min()
        ymax = y_data.max()
        axis = [xmin,xmax,ymin,ymax]
    else:
        xmin = axis[0]
        xmax = axis[1]
        ymin = axis[2]
        ymax = axis[3]

    # prepare data for jpdf plot
    X, Y      = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values    = np.vstack([x_data,y_data])
    kernel    = scipy.stats.gaussian_kde(values)
    Z         = np.reshape(kernel(positions).T, X.shape)
    
    
    plt.figure()
    plt.pcolor(X,Y,Z)
    plt.plot(x_data, y_data, 'k.', markersize=3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.axis([xmin, xmax,ymin,ymax])
    if not title == "":
        plt.title(title)
    #plt.set_ylim([ymin, ymax])
    cb = plt.colorbar()
    cb.set_label("probability density")
    plt.show()

def make_simulations_jpdf_plot(sim_results,
                               param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_all_sims[:,x_idx]
                y_data  = sim_results.np_all_sims[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot with performance constraints")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    
    
def make_cull_set_jpdf_plot(sim_results,
                            param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_pareto_set_cull[:,x_idx]
                y_data  = sim_results.np_pareto_set_cull[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot with performance constraints")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    
                    
def make_pareto_set_jpdf_plot(sim_results,
                              param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_pareto_set[:,x_idx]
                y_data  = sim_results.np_pareto_set[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    

def make_all_jpdf_plot(sim_results, param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_data  = sim_results.np_pareto_set[:,sim_results.all_names.index(x_label)]
                y_data  = sim_results.np_pareto_set[:,sim_results.all_names.index(y_label)]
                
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")
                
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                
                x_data  = sim_results.np_pareto_set_cull[:,x_idx]
                y_data  = sim_results.np_pareto_set_cull[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    pyflamestk.pareto.make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")

def resample_from_kernel_density(sim_results, 
                                 param_list, 
                                 param_list_not_free,
                                 fname_param = 'params.dat',
                                 n_resamples = 10000):


    # create free parameter list
    print("resampling from kernel density estimate...")
    print("\t output to file: {}".format(fname_param))
    print("\t number of resamples: {}".format(n_resamples))
    print("\t creating free parameter list...")
    free_param_list = param_list.copy()
    for param in param_list_not_free:
        free_param_list.remove(param)

    # print param list and whether or not the parameter is free        
    for param in param_list:
        if param in free_param_list:
            print("\t\t {} FREE".format(param))
        else:
            print("\t\t {} NOT FREE".format(param))

    # select only the parameters we are interested in
    param_col_index = []    
    for param_name in free_param_list:
        p_idx = sim_results.all_names.index(param_name)
        param_col_index.append(p_idx)


    #select only the parameters and do a KDE analysis
    print("\t do KDE analysis on culled parameter set...")
    cull_parameters = sim_results.np_pareto_set_cull[:,param_col_index]
    kernel = scipy.stats.gaussian_kde(cull_parameters.transpose())

    #resample from the kernel
    print("\t generation samples (may take a while)...")
    sim_results.param_sample = kernel.resample(size=n_resamples)

    #write file out
    param_list_out = []
    for idx, row in enumerate(sim_results.param_sample.T):
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
                    param_value = 0.
                elif param == 'p_MgMg_c':
                    param_value = 0.
                elif param == 'p_MgMg_rho':
                    param_value = 0.5
                elif param == 'p_MgO_c':
                    param_value = 0.
                else:
                    pass
            row_params.append(param_value)
        param_list_out.append(row_params)
    print("\t number of parameter sets: {}".format(len(param_list_out)))
    print("\t writing to file...")


    # print header
    str_out = "idx "
    for param in param_list:
      str_out += "{} ".format(param)
    
    # print params
    for idx,row in enumerate(param_list_out):
      str_out += "\n"
      str_out += "{} ".format(idx)
      for param_value in row:
        str_out += "{} ".format(param_value)
    
    
    f = open(fname_param, mode='w')
    f.write(str_out)  
    f.close()

def make_2d_pareto_plots(sim_results,
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
                                                   label='dominated',s=1))
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
    
                plt.axis([0,
                          np.percentile(sim_results.np_pareto_set[:,x_idx],90),
                          0,
                          np.percentile(sim_results.np_pareto_set[:,y_idx],90)])
                plt.legend(handles=plt_handles)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.show()

        
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],"l:s:")
    except getopt.GetoptError:
        print("pareto.py -l file1 -s file2 -l file3 ...")
        sys.exit(2)

    raw_vectors=[]
    for opt, arg in opts:
        if opt == '-l':
            raw_vectors.append(readfile(arg,1.0))
        elif opt == '-s':
            raw_vectors.append(readfile(arg,-1.0))

    dataset = create_dataset(raw_vectors)
    nondominated_sort(dataset)
    for s in dataset:
        print(s)
