import copy
import os.path
import numpy as np

class PyPosmatEngine:
    
  def __init__(self, 
               fname_config_pyposmat  = "pyposmat.config",
               fname_config_potential = "pyposmat.potential",
               fname_config_qoi       = "pyposmat.qoi",
               is_read = True):
                   
    self.supported_qoi = ['sp_ecoh','sp_pressure']
    print("pyflamestk.pyposmat version 0.1")
    self.read_configuration_files(fname_config_pyposmat=fname_config_pyposmat,
                                  fname_config_potential=fname_config_potential,
                                  fname_config_qoi=fname_config_qoi)
    self.check_lammps_configuration()
    self.check_structure_database()
    self.check_lammps_simulation_templates()
    self.check_external_software()
    self.check_quantities_of_interest()
    self.check_potential_parameters()

  def evaluate_parameter_set(self, potential_parameters):
      pass

  def write_lammps_potential_file(self):
      pass

  def create_lammps_simulation(self):
      pass

  def check_potential_parameters(self):
      
      self.potential_parameter_list = self.config_potential.get_param_list()
      self.potential_parameters = {}
      for param in self.potential_parameter_list:
          self.potential_parameters[param] = self.config_potential.get_param_value(param)
      
  def check_quantities_of_interest(self):
      print("\nchecking quantities of interest..")
      structures = []
      for qoi_key in self.config_qoi.qoi_info['qoi'].keys():
          new_structures = self.config_qoi.qoi_info['qoi'][qoi_key]['structure']
          for structure in new_structures:
              if structure not in structures:
                structures.append(structure)
                
      # confirm that structures requested in the quantities of interest
      # are in the fitting database
      for structure in structures:
          if structure in self.structure_db.keys():
              msg_out = "qoi_structure({}) in structure_db..."
              msg_out = msg_out.format(structure)
              print(msg_out)
          else:
              err_out = "qoi_structure({}) not in in structe_db."
              err_out = err_out.format(structure)
              raise ValueError(err_out)
              
      variables = []
      for qoi_key in self.config_qoi.qoi_info['qoi'].keys():
          new_variable = self.config_qoi.qoi_info['qoi'][qoi_key]['variable']
          if new_variable not in variables: 
              variables.append(new_variable)
              
              
      # confirm that the variables request in the quantities of interest
      # can be calculated with the known simulation templates
      for variable in variables:
          if variable in self.supported_qoi:
              msg_out = "qoi_variable({}) is in qoi_db"
              msg_out = msg_out.format(variable)
              print(msg_out)
          else:
              err_out = "qoi_variable({}) is not in qoi_db"
              err_out = err_out.format(variable)
              raise ValueError(err_out)

  def run(self):
      pass
  
  def check_external_software(self):
      pass
                                  
  def read_configuration_files(self,
                               fname_config_pyposmat, 
                               fname_config_potential,
                               fname_config_qoi):
      self.config_pyposmat  = ConfigFile(fname_config=fname_config_pyposmat, is_read = True)
      self.config_potential = PotentialConfigFile(fname_config=fname_config_potential, is_read = True)
      self.config_qoi       = QoiConfigFile(fname_config=fname_config_qoi)
      
  def check_structure_database(self):
      out_msg = "\nchecking if all LAMMPS structures templates exist..."
      print(out_msg)
      self.dir_structure_db = self.config_pyposmat.config_dict['dir_structure_db']
      self.structure_db = {}
      if not os.path.isdir(self.dir_structure_db):
          err_msg = "dir_structure_db({}) does not exist"
          err_msg = err_msg.format(self.dir_structure_db)
          raise IOError(err_msg)
      else:
          out_msg = "found dir_structure_db({})"
          out_msg = out_msg.format(self.dir_structure_db)
          print(out_msg)
          
      # get lammps structures
      for key in self.config_pyposmat.config_dict['structures'].keys():
          name   = key
          fname  = self.config_pyposmat.config_dict['structures'][key]['name']
          kind   = self.config_pyposmat.config_dict['structures'][key]['type']
          if (kind == 'lmmps'):
              if (os.path.isfile(os.path.join(self.dir_structure_db,fname))):
                  out_msg = "structure_db adding {} -> {}"
                  out_msg = out_msg.format(name,fname)
                  print(out_msg)
                  # if file exists, add to structure database
                  self.structure_db[name] = fname
              else:
                  # if file does not exist, do not add to structure database
                  err_msg = "file({}) for structure({}) does not exist"
                  err_msg = err_msg.format(fname,name)
                  raise IOError(err_msg)
          else:
              # if structure type is not supported, raise an error
              # structure types suppored = 'lmmps'
              err_msg = "structure_type({}) is not supported"
              err_msg = err_msg.format(kind)
              raise ValueError(err_msg)
  
  def check_lammps_configuration(self):
      # check lmps_bin
      # check lmps_exe_script
      pass
  
  def check_lammps_simulation_templates(self):
      out_msg = "\nchecking if all LAMMPS simulations templates exist..."
      print(out_msg)
      self.dir_lammps_sim_db = self.config_pyposmat.config_dict['dir_lammps_sim_db']
      self.lammps_sim_db = {}
      if not os.path.isdir(self.dir_lammps_sim_db):
          err_msg = "dir_lammps_sim_db({}) does not exist"
          err_msg = err_msg.format(self.dir_lammps_sim_db)
          raise IOError(err_msg)
      else:
          out_msg = "found dir_lammps_sim_db({})"
          out_msg = out_msg.format(self.dir_lammps_sim_db)
          print(out_msg)
          
      # get lammps simulation templates
      for key in self.config_pyposmat.config_dict['lmps_sim_type'].keys():
          name   = key
          fname  = self.config_pyposmat.config_dict['lmps_sim_type'][key]
          if (os.path.isdir(os.path.join(self.dir_lammps_sim_db,fname))):
              out_msg = "lammps_sim_db adding {} -> {}"
              out_msg = out_msg.format(name,fname)
              print(out_msg)
              # if file exists, add to structure database
              self.structure_db[name] = fname
          else:
              # if file does not exist, do not add to structure database
              err_msg = "directory({}) for sim_template({}) does not exist"
              err_msg = err_msg.format(fname,name)
              raise IOError(err_msg)

class MonteCarloParameterSampler(PyPosmatEngine):
    
    
  def run(self,n_simulations = 10):
    param_set = []  # initialize
    n_parameters = len(self.potential_parameter_list)
    for i_simulation in range(n_simulations):

        new_param_set = np.zeros(shape=[1,n_parameters])
        
        i = 0
        for param in self.potential_parameter_list:
            a = self.potential_parameters[param][1]
            b = self.potential_parameters[param][2]
            new_param_set[0,i] = np.random.uniform(low = a, high = b)
            i += 1
        self.evaluate_parameter_set(new_param_set[0,:].tolist())
        param_set.append(new_param_set[0,:].tolist())
    np_param_set = np.array(param_set)

  def evaluate_parameter_set(self, parameter_set):
      print(parameter_set)

class ConfigFile:
  def __init__(self, fname_config = "pyposmat.config",is_read = True):
    self.fname_config = fname_config
    if is_read == True:
      self.read()
      
  def get_pyposmat_info(self):
      self.pyposmat_type = self.config_dict['pyposmat_type']
      self.pyposmat_log  = self.config_dict['pyposmat_log_fname']
      self.pyposmat_out  = self.config_dict['pyposmat_out_fname']
      return {'pyposmat_type': self.pyposmat_type,
              'pyposmat_log' : self.pyposmat_log,
              'pyposmat_out' : self.pyposmat_out}

  def get_number_of_objective_functions(self):
    n_qoi         = len(self.config_dict['qoi'])
    n_qoi_weights = len(self.config_dict['qoi_weights'])
    n_qoi_targets = len(self.config_dict['qoi_target'])
    
    if (n_qoi == n_qoi_targets) and (n_qoi == n_qoi_targets):
        # check to see if qoi, weights, and values are equal
        return n_qoi
    else:
        err_msg = 'n_qoi({}), n_qoi_weights({}), n_qoi_targets({}) are different.'
        err_msg.format(n_qoi, n_qoi_weights, n_qoi_targets)
        raise ValueError(err_msg)
      
  def getNumberOfObjectiveFunctions(self):
    n_qoi         = len(self.config_dict['qoi'])
    n_qoi_weights = len(self.config_dict['qoi_weights'])
    n_qoi_targets = len(self.config_dict['qoi_target'])
    return n_qoi

  def getDakotaKwargs(self):
    dakota_kwargs = {'dakota_bin'      :self.config_dict["dakota_bin"],
                     'dakota_fname_in' :self.config_dict["dakota_fname_in"],
                     'dakota_fname_out':self.config_dict["dakota_fname_out"],
                     'dakota_fname_err':self.config_dict["dakota_fname_err"]}
    return dakota_kwargs

  def read(self):
    self.config_dict = read_pyposmat_config_file(self.fname_config)
    return copy.deepcopy(self.config_dict)
    
  
def read_pyposmat_config_file(fname_config):
  config_dict = {}
  config_dict['structures'] = {}
  config_dict['lmps_sim_type'] = {}
  config_dict['qoi_list'] = []
  config_dict['qoi'] = {}
  config_dict['qoi_weights'] = {}
  config_dict['qoi_targets'] = {}
  config_dict['qoi_normfactor'] = {}
  file = open(fname_config,'r')
  for idx, line in enumerate(file.readlines()):
    if line.strip().startswith('#'):
      # skip if the line is a comment
      pass
    elif line.strip() == '':
      # skip if the line is empty or only contains white space
      pass
    else:
      config_info = line.split('=')
      if len(config_info) > 0:
        config_var = config_info[0].strip()
        config_param = config_info[1].strip().strip('"').strip()
        if config_var == 'structure':
          # structure = structure_name, structure_fname, structure_type
          structure_name  = config_param.split(',')[0].strip()
          structure_fname = config_param.split(',')[1].strip()
          structure_type  = config_param.split(',')[2].strip()
          config_dict['structures'][structure_name] = {}
          config_dict['structures'][structure_name]['name'] = structure_fname
          config_dict['structures'][structure_name]['type'] = structure_type
        elif config_var == "lmps_sim_type":
          sim_type     = config_param.split(',')[0].strip()
          sim_location = config_param.split(',')[1].strip()
          config_dict['lmps_sim_type'][sim_type] = sim_location
        elif config_var == "qoi_list":
          config_dict["qoi_list"] = [word.strip() for word in config_param.split(',')]
        elif config_var == "qoi":
          args = config_param.split(',')
          n_structures = len(args) - 2
          qoi_name = args[0].strip()
          qoi_type = args[1].strip()
          config_dict['qoi'][qoi_name] = {}
          config_dict['qoi'][qoi_name]['type'] = qoi_type
          config_dict['qoi'][qoi_name]['structures'] = [args[i].strip() for i in range(2,2 + n_structures)]
        elif config_var == "qoi_weights":
          qoi_name   = config_param.split(',')[0].strip()
          qoi_weight = float(config_param.split(',')[1].strip())
          config_dict['qoi_weights'][qoi_name] = qoi_weight
        elif config_var == "qoi_target":
          qoi_name       = config_param.split(',')[0].strip()
          qoi_target_val = float(config_param.split(',')[1].strip())
          config_dict['qoi_targets'][qoi_name] = qoi_target_val
        elif config_var == "qoi_normfactor":
          qoi_name       = config_param.split(',')[0].strip()
          qoi_normf      = float(config_param.split(',')[1].strip())
          config_dict['qoi_normfactor'][qoi_name] = qoi_normf
        else:
          config_dict[config_var] = config_param
  return config_dict

class PotentialConfigFile:
  def __init__(self, fname_config = "pyposmat.potential",is_read = True):
    self.fname_config = fname_config
    self.config_dict  = {}
    if is_read == True:
      self.read()
      
  def read(self):
    self.elements = {}
    self.pair     = {}
    file = open(self.fname_config)
    for idx, line in enumerate(file.readlines()):
        if line.strip().startswith('#'):
            pass
        elif line.strip() == '':
            pass
        else:
            keyword = line.split('=')[0].strip()
            params  = line.split('=')[1].strip().split()
            if keyword == 'potential_elements':
                for param in params:
                    self.elements[param.strip()] = {}
            elif keyword == 'potential_charge':
                element  = params[0].strip()
                chrg_0    = float(params[1])
                chrg_min  = float(params[2])
                chrg_max  = float(params[3])
                if chrg_min > chrg_max:
                    raise ValueError("chrg_min < chrg_max failed. chrg_min = {}. chrg_max = {}".format(chrg_min,chrg_max))
                self.elements[element]['charge'] = [chrg_0, chrg_min, chrg_max]
            elif keyword == 'potential_pair_type':
                element1  = params[0].strip()
                element2  = params[1].strip()
                pair_type = params[2].strip()
                pair_str  = "{}{}".format(element1,element2)
                self.pair[pair_str] = {}
                self.pair[pair_str]['type'] = pair_type
                self.pair[pair_str]['param'] = {}
            elif keyword == 'potential_pair_param':
                element1   = params[0].strip()
                element2   = params[1].strip()
                param_name = params[2].strip()
                param_0    = float(params[3])
                param_min  = float(params[4])
                param_max  = float(params[5])
                pair_str   = "{}{}".format(element1,element2)
                self.pair[pair_str]['param'][param_name] = [param_0,param_min,param_max]
                              
  def get_param_list(self):
    self.param_list = []
    for element in self.elements:
        self.param_list.append('chrg_{}'.format(element,self.elements[element]['charge']))
    for pair in self.pair:
        for param in self.pair[pair]['param']:
            self.param_list.append('p_{}_{}'.format(pair,param))
    return self.param_list

  def get_param_value(self,param_name):
      pot_type = param_name.split('_')[0]
      if pot_type == 'p':
          pot_id   = param_name.split('_')[1]
          pot_param = param_name.split('_')[2]
          return self.pair[pot_id]['param'][pot_param]
      elif pot_type == 'chrg':
          pot_id   = param_name.split('_')[1]
          return self.elements[pot_id]['charge']
          
class QoiConfigFile:
  def __init__(self, fname_config = "pyposmat.potential",is_read = True):
    self.fname_config = fname_config
    if is_read == True:
      self.read()
      
  def read(self):
    self.qoi_info = {}
    self.qoi_info['qoi'] = {}
    file = open(self.fname_config)
    for idx, line in enumerate(file.readlines()):
        if line.strip().startswith('#'):
            pass
        elif line.strip() == '':
            pass
        else:
            keyword = line.split('=')[0].strip()
            params  = line.split('=')[1].strip()
            if keyword == 'qoi_list':
                qoi_list = params.split(',')
                self.qoi_info[keyword] = qoi_list
            elif keyword == 'define_qoi':
                qoi_list      = params.split(',')
                
                qoi_name      = qoi_list[0]
                qoi_variable  = qoi_list[1].strip()
                qoi_structures = [qoi_list[i].strip() for i in range(2,len(qoi_list))]
                
                self.qoi_info['qoi'][qoi_name] = {}
                self.qoi_info['qoi'][qoi_name]['variable']  = qoi_variable
                self.qoi_info['qoi'][qoi_name]['structure'] = qoi_structures                
            elif keyword == 'qoi_target':
                qoi_list     = params.split(',')
                
                qoi_name     = qoi_list[0].strip()
                qoi_target   = float(qoi_list[1].strip())
                qoi_weight   = float(qoi_list[2].strip())
                qoi_norm     = float(qoi_list[3].strip())
                
                self.qoi_info['qoi'][qoi_name]['target'] = qoi_target
                self.qoi_info['qoi'][qoi_name]['weight'] = qoi_weight
                self.qoi_info['qoi'][qoi_name]['norm']   = qoi_norm
                
                print(self.qoi_info['qoi'][qoi_name])
            else:
                msg_err = "unknown qoi_keyword({})"
                msg_err = msg_err.format(keyword)
                raise ValueError(msg_err)