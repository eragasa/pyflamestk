import sys
import re
import yaml

class DakotaInterface:

        def __init__(self,
                     filename_param,
                     filename_results):
            self.filename_parameters = filename_param
            self.filename_results    = filename_results
            
            self.yaml_out = "yaml_out.dat"            
            self.parameter_dict = {}
            self.eval_id = 0
            self.num_functions = 0
            self.num_variables = 0
            self.num_deriv_variables = 0
            self.num_analysis_components = 0
            self.variables = {}
            self.active_state_vectors = {}
            self.deriv_values_vectors = {}
            self.active_components = {}

            self._readParametersFile()
            self._parseParametersDictionary()
            self._checkParameters()
            self.sendToApplication()
            self._sendResultsToDakota()

        def sendToApplication(self):
            self._writeDakotaOut()
            # execute the rosenbrock analysis as a separate Python module
            from rosenbrock import RosenbrockSimulation
            rs = RosenbrockSimulation(self.yaml_out_str, "results.out.tmp")


        def _readParametersFile(self):

            # setup regular expressions for parameter/label matching
            e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
            f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
            i = '-?\\d+'                                         # integer
            value = e+'|'+f+'|'+i                                # numeric field
            tag = '\\w+(?::\\w+)*'                               # text tag field
            # regular expression for aprepro parameters format
            aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
            # regular expression for standard parameters format
            standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')

            file_parameters = open(self.filename_parameters)
            self.parameter_dict = {}
            for line in file_parameters:
                m = aprepro_regex.match(line)
                if m:
                    self.parameter_dict[m.group(1)] = m.group(2)
                else:
                    m = standard_regex.match(line)
                    self.parameter_dict[m.group(2)] = m.group(1)
            file_parameters.close()
            
        def _parseParametersDictionary(self):
            
            for key in self.parameter_dict.keys():
                if key == "eval_id":
                    self.eval_id = int(self.parameter_dict['eval_id'])
                # procesing functions
                elif key == "functions":
                    self.num_functions = int(self.parameter_dict[key])
                elif key == "DAKOTA_FNS":
                    self.num_functions = int(self.parameter_dict[key])
                elif key.startswith("ASV"):
                    # ASV = 1 = dakota requesting response_function value
                    # ASV = 2 = dakota requesting response function gradient
                    # ASV = 4 = dakota requesting response function hessian
                    self.active_state_vectors[key] = int(self.parameter_dict[key])
                # processing derivative values vector
                elif key == "derivative_variables":
                    self.num_deriv_variables = int(self.parameter_dict[key])
                elif key.startswith("DVV"):
                    self.deriv_values_vectors[key] = int(self.parameter_dict[key])
                # analysis components
                elif key == "analysis_components":
                    self.num_analysis_components = int(self.parameter_dict[key])
                elif key.startswith("AC"):
                    self.active_components = self.parameter_dict[key]
                # variables
                elif key == "variables":
                    self.num_variables = int(self.parameter_dict['variables'])
                elif key == "DAKOTA_VARS":
                    self.num_variables = int(self.parameter_dict['variables'])
                else:
                    # everything must be a variable
                    self.variables[key] = float(self.parameter_dict[key])

        def _checkParameters(self):
            # crude error checking; handle both standard and aprepro cases
            self.num_variables = 0
            if ('variables' in self.parameter_dict):
                self.num_variables = int(self.parameter_dict['variables'])
            elif ('DAKOTA_VARS' in self.parameter_dict):
                self.num_variables = int(self.parameter_dict['DAKOTA_VARS'])
            
            self.num_functions = 0
            if ('functions' in self.parameter_dict):
                self.num_functions = int(self.parameter_dict['functions'])
            elif ('DAKOTA_FNS' in self.parameter_dict):
                self.num_functions = int(self.parameter_dict['DAKOTA_FNS'])
            
            if (self.num_variables != 2 or self.num_functions != 1):
                strErr = "Rosenbrock requires 2 variables and 1 function;\nfound {} variables and {} functions."
                print(strErr.format(str(self.num_variables), str(self.num_functions)))
                sys.exit(1)

        def _writeDakotaOut(self):
            data = dict(
                eval_id       = self.eval_id,
                num_variables = self.num_variables,
                variables     = self.variables,
                num_functions = self.num_functions,
                asv           = self.active_state_vectors,
                dvv           = self.deriv_values_vectors,
                ac            = self.active_components
                )
            self.yaml_out_str = yaml.dump(data, default_flow_style=False)
            file = open(self.yaml_out,'w')
            file.write(self.yaml_out_str)
            file.close()


        def _sendResultsToDakota(self):
            # move the temporary results file to the one DAKOTA expects
            import shutil
            shutil.move('results.out.tmp', self.filename_results)
            #os.system('mv results.out.tmp ' + sys.argv[2])
