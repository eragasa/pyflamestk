import copy
import os
import os.path
import shutil
import subprocess
import numpy as np

import pyflamestk.lammps
class PyPosmatError(Exception):
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)

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
                                  
    # list of class member variables
    self.qoi = {}
    self.lammps_sims = {}        #lammps_sims[sim_name]['structure','sim_type']
    self.potential_parameter_list = []
    self.potential_parameters = {}
    self.dir_lammps_sim_db = None
    
    self.check_lammps_configuration()
    self.check_structure_database()
    self.check_lammps_simulation_templates()
    self.check_external_software()
    self.check_quantities_of_interest()
    self.determine_lammps_simulations()
    self.check_potential_parameters()

    #self.write_lammps_potential_file()

  def evaluate_parameter_set(self, potential_parameters):
      potential_type = 'buckingham'        # TODO: GET RID OF MAGIC VARIABLE

      fname_potential = None  # initialize, default parameterization requires None
      fname_eam       = None  # initialize, default parameterization requires None

      fname_potential = 'potential.mod'

      if(potential_type == 'buckingham'):
          fname_potential = 'potential.mod'
          self.write_lammps_buckingham_potential_file(potential_params=potential_parameters, 
                                                      fname_out = fname_potential)
      elif(potential_type == 'eam'):
          fname_potential = "potential.mod"
          fname_eam       = "eam.alloy"
          self.write_eam_potential_file()
          self.write_eam_setfl_file()
      else:
          msg_err = "unknown potential_type code({})"
          msg_err = msg_err.format(potential_type)
          raise ValueError(msg_err)          
          
      #self.create_lammps_simulations()
      for simulation in self.lammps_sims:
          structure_name = simulation.split('.')[0]
          sim_code = simulation.split('.')[1]
          
          self.lammps_sims[simulation]['structure'] = structure_name
          
          if (sim_code == 'sp'):
              self.lammps_sims[simulation]['sim_type']  = 'single_point'
          else:
              msg_err = "unknown sim_type code({})"
              msg_err = msg_err.format(sim_code)
              raise ValueError(msg_err)

      for simulation in self.lammps_sims:
          name_structure  = self.lammps_sims[simulation]['structure']
          fname_structure = self.structure_db[name_structure]
          self.create_lammps_simulation(sim_name=simulation,
                                        fname_structure=fname_structure,
                                        fname_sim_template=self.lammps_sims[simulation]['sim_type'],
                                        fname_potential = fname_potential,
                                        fname_eam = fname_eam)

  def run_all_lammps_simulations(self):
    #TODO: replace this line with information in the config file.
    fname_script_run_lammps = "./runsimulation.sh"
    for dir in self.lammps_sims:
      print("  ",dir)
      p = subprocess.Popen(fname_script_run_lammps, shell=True, cwd=dir)

    sims_finished  = False
    sim_has_error = False

    while not(sims_finished) and not(sim_has_error):
      sims_finished  = self.check_all_lammps_simulations_done(sim_directories = self.lammps_sims)
      sim_has_error  = self.check_all_lammps_simulations_for_errors(sim_directories = self.lammps_sims) 

    if(sim_has_error):
      print("simulation has failed")      

  def check_all_lammps_simulations_done(self, sim_directories):
    return pyflamestk.lammps.checkAllLammpsSimulationsDone(sim_directories = sim_directories)

  def check_all_lammps_simulations_for_errors(self, sim_directories):
    return pyflamestk.lammps.checkAllLammpsSimulationsForErrors(sim_directories = sim_directories)

  def determine_lammps_simulations(self):
      
      ''' This method determines which lammps simulations are required from
      the QOI simulation list.'''
            
      self.lammps_sims = {}   # initialization
      print("determining lammps simulations...")
      for qoi in self.qoi_list:
          if self.qoi[qoi]['variable'] in ['sp_ecoh', 'sp_pressure']:
              structure = self.qoi[qoi]['structure'][0]
              lammps_sim_template = "sp"
              lammps_sim_name = "{}.{}"
              lammps_sim_name = lammps_sim_name.format(structure,lammps_sim_template)
              if lammps_sim_name not in self.lammps_sims.keys():
                  self.lammps_sims[lammps_sim_name] = {}
                  print("  ",lammps_sim_name, "created.")
          else:
              # qoi is not supported
              msg_err = "qoi({}) does not support type({})..."
              msg_err = msg_err.format(qoi,self.qoi[qoi]['variable'])
              raise ValueError(msg_err)

  def get_mass(self,element):
      if element == 'Mg':
          return 24.305
      elif element == "O":
          return 15.999
      else:
          raise ValueError("element {} not in database".format(element))
          
  def get_name(self,element):
      if element == "Mg":
          return 'magnesium'
      elif element == "O":
          return 'oxygen'
      else:
          raise ValueError('element {} not in database'.format(element))
          
  def write_lammps_buckingham_potential_file(self, potential_params, fname_out):
      elements = self.potential_elements.keys()
      charge = {'Mg':+2.00,'O':-2.00}
      str_pot = ""
      for i, element in enumerate(elements):
          str_pot += "mass {} {}\n".format(i+1,self.get_mass(element))
      str_pot += "\n"
      for i, element in enumerate(elements):
          str_pot += "group {} type {}\n".format(self.get_name(element),i+1)
      str_pot += "\n"
      for i, element in enumerate(elements):
          str_pot += "set group {} charge {}\n".format(self.get_name(element),
                                                       charge[element])
      str_pot += "\n"
      str_pot += "variable R_cut equal 10.0\n"
      str_pot += "\n"
      str_pot += "pair_style buck/coul/long ${R_cut}\n"
      
      for i, element_i in enumerate(elements):
          for j, element_j in enumerate(elements):
              if i <= j:
                  try:
                      p_rho = potential_params['p_{}{}_rho'.format(element_i,element_j)]
                      p_a   = potential_params['p_{}{}_a'.format(element_i,element_j)]
                      p_c   = potential_params['p_{}{}_c'.format(element_i,element_j)]
                  except KeyError:
                      
                      p_rho = 0.0
                      p_a = 0.5
                      p_c = 0.0
                  str_pot += "pair_coeff {} {} {} {} {} {}\n".format(i+1,j+1,p_rho,p_a,p_c,'${R_cut}')
      str_pot += "\n"
      str_pot += "kspace_style pppm 1.0e-5\n"
      str_pot += "\n"
      str_pot += "neighbor 1.0 bin\n"
      str_pot += "neigh_modify every 1 delay 0 check yes\n"
      f = open(fname_out,'w')
      f.write(str_pot)
      f.close()
      return str_pot

  def write_lammps_potential_file(self): 
      self.write_lammps_buckingham_potential_file()      

  def create_lammps_simulations(self):
      for sim_name in self.lammps_sims:
          print(sim_name)
          
  def create_lammps_simulation(self,
                               sim_name,
                               fname_structure,
                               fname_sim_template,
                               fname_potential = None,
                               fname_eam = None):
      potential_type = 'buckingham'
      dir_name = sim_name
      
      # remove the directory if it exists
      if os.path.exists(dir_name):
          shutil.rmtree(dir_name)

      # copy simulation template
      src_dir = os.path.join(self.dir_lammps_sim_db, fname_sim_template)
      dst_dir = os.path.join(os.path.curdir,dir_name)
      shutil.copytree(src=src_dir,dst=dst_dir)
      
      # copy lammps structure file
      src_fname = os.path.join(self.dir_structure_db, fname_structure)
      dst_fname = os.path.join(os.path.curdir,dir_name,'lammps.structure')
      shutil.copyfile(src= src_fname, dst=dst_fname)
      
      if not(fname_potential == None):
        src_fname = os.path.join(fname_potential)
        dst_fname = os.path.join(os.path.curdir,dir_name,'potential.mod')
        shutil.copyfile(src_fname,dst_fname)
      else:
        print('not copying {}'.format(fname_potential))
      
      if not(fname_eam == None):
          shutil.copyfile(fname_eam,"{}/eam.alloy".format(dir_name))

      #modify input file
      if potential_type == 'buckingham':
        # get the name of the input file
        path_sim = os.path.join(os.path.curdir,dir_name)
        fname_list = [f for f in os.listdir(path_sim) if os.path.isfile(os.path.join(path_sim,f))]
        fname_input_file = ""
        for fname in fname_list:
          if fname.startswith('in.'):
            fname_input_file = fname
        # replace atom type
        # step 1: write new input file with modifications
        # step 2: delete old input file
        # step 3: move new input file into old input file location
        old_fname = os.path.join(path_sim,fname_input_file)
        new_fname = os.path.join(path_sim,"{}.tmp".format(fname_input_file))
        with open(new_fname,'w') as new_file:
          with open(old_fname,'r') as old_file:
            for line in old_file:
              new_file.write(line.replace('atom_style atomic', 'atom_style charge'))
        os.remove(old_fname)             # remove old input file
        shutil.move(new_fname,old_fname)  # copy new input file to old input file

  def check_potential_parameters(self):

      self.potential_elements = self.config_potential.elements      
      self.potential_parameter_list = self.config_potential.get_param_list()
      self.potential_parameter_types = {}
      for pair in self.config_potential.pair.keys():
          self.potential_parameter_types[pair] = self.config_potential.pair[pair]['type']
      self.potential_parameters = {}
      for param in self.potential_parameter_list:
          self.potential_parameters[param] = self.config_potential.get_param_value(param)
      self.potential_types = self.config_potential.pair
      
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
              err_out = "qoi_structure({}) not in in structure_db."
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

      # build qoi lists      
      self.qoi_list = self.config_qoi.qoi_info['qoi'].keys()
      self.qoi = {}          
      for qoi_key in self.qoi_list:
          self.qoi[qoi_key] = {}
          self.qoi[qoi_key]['structure']   = self.config_qoi.qoi_info['qoi'][qoi_key]['structure']
          self.qoi[qoi_key]['weight']      = self.config_qoi.qoi_info['qoi'][qoi_key]['weight']
          self.qoi[qoi_key]['target']      = self.config_qoi.qoi_info['qoi'][qoi_key]['target']
          self.qoi[qoi_key]['norm']        = self.config_qoi.qoi_info['qoi'][qoi_key]['norm']
          self.qoi[qoi_key]['variable']    = self.config_qoi.qoi_info['qoi'][qoi_key]['variable']
          print(qoi_key)
          print("\t",self.qoi[qoi_key]['structure'])
          print("\t",self.qoi[qoi_key]['weight'])
          print("\t",self.qoi[qoi_key]['target'])
          print("\t",self.qoi[qoi_key]['norm'])
          print("\t",self.qoi[qoi_key]['variable'])
         
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

  def calculate_qoi(self):
    self.get_properties_from_simulations()
    self.get_qoi_estimates()

    self.abs_error = {}
    self.square_error = {}
    self.square_error_normalized = {}

    for key in self.qoi_list:
      self.abs_error[key] = abs(self.qoi[key]['est'] - self.qoi[key]['target'])
      self.square_error[key] = self.abs_error[key]**2
      self.square_error_normalized[key] = self.square_error[key]/self.qoi[key]['norm']

  def get_properties_from_simulations(self):
    # get properties from simulation
    self.properties = {}
    for qoi_key in self.qoi_list:
      variable_name = self.qoi[qoi_key]['variable']
      structures    = self.qoi[qoi_key]['structure']

      if variable_name == 'latt_a':
        simulation      = "{}.{}".format(structures[0],'E_coh')
        property_name   = "{}.{}".format(structures[0],'a_lat')
        if not(property_name in self.properties.keys()):
          self.properties[property_name] = pyflamestk.lammps.getLatticeParameter(os.path.join(simulation,'out.dat'))
      elif variable_name == 'E_coh':
        simulation      = "{}.{}".format(structures[0],'E_coh')
        property_name   = "{}.{}".format(structures[0],'E_coh')
        if not(property_name in self.properties.keys()):
          self.properties[property_name] = pyflamestk.lammps.getCohesiveEnergy(os.path.join(simulation,'out.dat'))
      elif variable_name == 'phase_difference':
        for structure in structures:
          simulation    = "{}.{}".format(structure,'E_coh')
          property_name = "{}.{}".format(structure,'E_coh')
          self.properties[property_name] = pyflamestk.lammps.getCohesiveEnergy(os.path.join(simulation,'out.dat'))
      elif variable_name in ['c11', 'c12', 'c44', 'B_modulus', 'G_modulus']:
        simulation      = "{}.{}".format(structures[0],'elastic')
        for var in ['c11','c12','c44']:
          property_name   = "{}.{}".format(structures[0],var)
          self.properties[property_name]   = pyflamestk.lammps.getElasticComponent(os.path.join(simulation,'out.dat'),var)
      elif variable_name == "E_f_defect":
        for structure in structures:
          simulation    = "{}.{}".format(structure,'E_coh')
          property_name = "{}.{}".format(structure,'E_coh')
          self.properties[property_name] = pyflamestk.lammps.getCohesiveEnergy(os.path.join(simulation,'out.dat'))
      elif variable_name == "sp_pressure":
        simulation      = "{}.{}".format(structures[0],"sp")
        property_name   = "{}.{}".format(structures[0],'pressure_total')
        self.properties[property_name]   = pyflamestk.lammps.getPressure(os.path.join(simulation,'out.dat'),type='total')
      elif variable_name == "sp_ecoh":
        simulation      = "{}.{}".format(structures[0],"sp")
        property_name   = "{}.{}".format(structures[0],'E_coh_press')
        self.properties[property_name]   = pyflamestk.lammps.getCohesiveEnergy(os.path.join(simulation,'out.dat'))
      else:
        print("pyposmat qoi not supported: {}".format(variable_name))

  def get_qoi_estimates(self):
    for qoi_key in self.qoi_list:
      variable_name = self.qoi[qoi_key]['variable']
      structures    = self.qoi[qoi_key]['structure']

      if variable_name == 'latt_a':
        self.qoi[qoi_key]['est'] =  self.properties["{}.{}".format(structures[0],'a_lat')]
      elif variable_name == 'E_coh':
        self.qoi[qoi_key]['est'] =  self.properties["{}.{}".format(structures[0],'E_coh')]
      elif variable_name == 'phase_difference':
        property_names   = ["{}.{}".format(structure, 'E_coh') for structure in structures]
        e_coh_0 = self.properties[property_names[0]] # hcp
        e_coh_1 = self.properties[property_names[1]] # fcc
        print(property_names[0],e_coh_0)
        print(property_names[1],e_coh_1)
        self.qoi[qoi_key]['est'] = e_coh_0 - e_coh_1
      elif variable_name in ['c11', 'c12', 'c44']:
        self.qoi[qoi_key]['est'] = self.properties["{}.{}".format(structures[0],variable_name)]
      elif variable_name in ['B_modulus']:
        c11 = self.properties["{}.c11".format(structures[0])]
        c12 = self.properties["{}.c12".format(structures[0])]
        self.qoi[qoi_key]['est'] = calculate_bulk_modulus(c11,c12)
      elif variable_name in ['G_modulus']:
        c11 = self.properties["{}.c11".format(structures[0])]
        c12 = self.properties["{}.c12".format(structures[0])]
        self.qoi[qoi_key]['est'] = calculate_shear_modulus(c11,c12)
      elif variable_name == "E_f_defect":
        n_atoms_ideal  = 1 
        n_atoms_defect = 2
        E_coh_ideal  = self.properties["{}.{}".format(structures[0],'E_coh')]
        E_coh_defect = self.properties["{}.{}".format(structures[1],'E_coh')]
        self.qoi[qoi_key]['est'] = n_atoms_defect * E_coh_defect - n_atoms_defect * E_coh_ideal 
      elif variable_name == "sp_pressure":
        property_name   = "{}.{}".format(structures[0],'pressure_total')
        self.qoi[qoi_key]['est'] = self.properties[property_name]
      elif variable_name == "sp_ecoh":
        property_name   = "{}.{}".format(structures[0],'E_coh_press')
        self.qoi[qoi_key]['est'] = self.properties[property_name]
      else:
        print("pyposmat qoi not supported: {}".format(variable_name))

class MonteCarloParameterSampler(PyPosmatEngine):
    
  def run(self,n_simulations = 10, fname_results = 'results.out'):

    #check to see if file exists and delete
    if os.path.exists(fname_results):
      if os.path.isfile(fname_results):
        os.remove(fname_results)
    
    f = open(fname_results,'w')
    
    # print header line
    str_out = ""
    str_out += "sim_id "
    for param in self.potential_parameter_list:
      str_out += "{} ".format(param)
    str_out += "| "
    for qoi_key in self.qoi_list:
      str_out += "{}_abserr ".format(qoi_key)
      str_out += "{}_sqerr ".format(qoi_key)
      str_out += "{}_nsqerr ".format(qoi_key)
    str_out += "\n"
    f.write(str_out)

    param_set = []  # initialize
    qoi_set   = []  # initialize
    for i_simulation in range(n_simulations):

        new_param_set = {}
        i = 0
        for param in self.potential_parameter_list:
            a = self.potential_parameters[param][1]
            b = self.potential_parameters[param][2]
            new_param_set[param] = np.random.uniform(low = a, high = b)
            i += 1

        # TODO: generalize this code
        new_param_set['chrg_O'] = - new_param_set['chrg_Mg']

        print('evaluating param_set_id:{}'.format(i_simulation))
        self.evaluate_parameter_set(new_param_set)
        self.run_all_lammps_simulations()
        self.calculate_qoi()

        new_qoi_set = []
        i = 0
        for qoi_key in self.qoi_list:
          new_qoi_set.append(self.abs_error[qoi_key])
          new_qoi_set.append(self.square_error[qoi_key])
          new_qoi_set.append(self.square_error_normalized[qoi_key])

        str_out = ""
        str_out = "{} ".format(i_simulation)
        for param in self.potential_parameter_list:
          str_out += "{} ".format(new_param_set[param])
        str_out += "| "
        for qoi_result in new_qoi_set:
          str_out += "{} ".format(qoi_result)
        str_out += "\n"
        f.write(str_out)
        param_set.append(new_param_set)
        qoi_set.append(new_qoi_set)

    f.close()
    str_out = ""
    for i_set, set in enumerate(param_set):
      str_out += "{} ".format(i_set)
      for param in self.potential_parameter_list:
        str_out += "{} ".format(set[param])
      str_out += "\n"

    #np_param_set = np.array(param_set)

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

    # class member variables
    self.fname_config = fname_config
    self.config_dict  = {}
    self.elements = {}
    self.pair = {}
    
    
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
                
            else:
                msg_err = "unknown qoi_keyword({})"
                msg_err = msg_err.format(keyword)
                raise ValueError(msg_err)
