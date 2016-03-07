class MonteCarloParameterSampler:
  def __init__(self, fname_config = "pyposmat.config", is_read = True):
      
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