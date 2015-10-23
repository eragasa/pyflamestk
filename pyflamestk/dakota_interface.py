import sys
import re
import yaml

class InputFile:
    def __init__(self,fname,fname_analysis_driver):
        self.fname_in  = fname
        self.fname_out = fname.split(".")[0] + ".out"
        self.fname_err = fname.split(".")[0] + ".err"
        self.fname_analysis_driver = ""
        self.isCheck = True
        self.fname_outfile = "" # filename overidden by command line option
        self.output = "normal"
        self.variables = {}
        self.cv_keys = []       # continuous variable keys
        self.n_cv = 0           # number of continuous variables
        
    def addVariable(self, descriptor, obj_variable):
        self.variables[descriptor] = obj_variable
        
    def write(self):
        file = open(self.name, "w")
        file.write(self.getEnvironmentBlock())
        file.write(self.getInterfaceBlock())
        file.write(self.getMethodBlock())
        file.write(self.getModelBlock())
        file.write(self.getVariablesBlock())
        file.close()
        
    def getEnvironmentBlock(self):
        if self.fname_out == "":
            self.fname_out = self.fname.split(".")[0] + ".out"
        if self.fname_err == "":
            self.fname_out == self.fname.split(".")[0] + ".err"

        strOut = "environment\n"
        if self.isCheck:
            strOut += "\tcheck\n"
        strOut += "\toutput_file {}\n".format(self.fname_out)
        strOut += "\terror_file {}\n".format(self.fname_err)
        
        return strOut
    
    def getInterfaceBlock(self):
        self.fname_params  = "params.in"
        self.fname_results = "results.out"
        strOut  = "interface\n"
        strOut += "\tanalysis_drivers \'{}\'\n".format(fname_exe)
        strOut += "\t\tfork\n"
        strOut += "\t\t\tparameters_file = \'{}\'\n".format(self.fname_params)
        strOut += "\t\t\tresults_file = \'{}\'\n".format(self.fname_results)
        strOut += "\t\t\tfile_tag\n"
        strOut += "\t\t\tfile_save\n"
        
    def getMethodBlock(self):
        strOut = "method\n"
        return strOut
        
    def getModelBlock(self):
        strOut = "model\n"
        return strOut
    
    def getContinuousVariablesSubblock(self):
        self.cv_keys = []
        self.n_cv = 0

        for key in self.variables.keys():
            if self.variables[key].variable_type == "continuous_design":
                self.cv_keys.append(key)

        self.n_cv = len(self.cv_keys)

        if self.n_cv == 0:
            return ""

        lower_bounds = []
        upper_bounds = []
        descriptors  = []
        initial_point = []
        scale_types = []
        scale_values = []
        for key in self.variables.keys():
            descriptors.append(self.variables[key].descriptor)
            lower_bounds.append(self.variables[key].lower_bound)
            upper_bounds.apppend(self.variables[key].upper_bound)
            initial_point.append(self.variables[key].initial_point)
            scale_types.append(self.variables[key].scale_types)
            scale_values.append(self.variables[key].scale_values)
            
        str_descriptors  = "descriptors"
        for desc in descriptors:
            str_descriptors += " \'{}\'".format(desc)

        str_lower_bounds = "lower_bounds"
        for val in lower_bounds:
            str_lower_bounds += " {:10.6f}".format(val)

        str_upper_bounds = "upper_bounds"
        for val in upper_bounds:
            str_upper_bounds += " {:10.6f}".format(val)

        str_initial_point = "initial_point"
        for val in initial_point:
            str_initial_point += " {:10.6f}".format(val)
        
        str_scale_types = "scale_types ="
        for val in scale_types:
            str_scale_types += " \'{}\'".format(val)
            
        str_scale_values = "scales ="
        for vale in scale_values:
            str_scale_values += " \'{}\'".format(val)
            
        strOut  = "\tcontinuous_design = {}\n".format(self.n_cvs)
        strOut += "\t\t{}\n".format(str_descriptors)
        strOut += "\t\t{}\n".format(str_lower_bounds)
        strOut += "\t\t{}\n".format(str_upper_bounds)
        strOut += "\t\t{}\n".format(str_initial_point)
        strOut += "\t\t{}\n".format(str_scale_types)
        strOut += "\t\t{}\n".format(str_scale_values)
        return strOut
    
    def getVariablesBlock(self):
        strOut  = "variables\n"
        strOut += self.getContinuousVariablesSubblock()

class SensitivityAnalysisInputFile(InputFile):
    def __init__(self,fname):
        InputFile.__init__(self,fname)


class Variable:
    def __init__(self, descriptor, var_type):
        self.descriptor
        self.variable_type = var_type
        

class ContinuousDesignVariable(Variable):
    def __init__(self,
                 descriptor,
                 initial_point,
                 lower_bound,
                 upper_bound,
                 scale_type = "none",
                 scale_value = 0):
        Variable.__init__(self,
                          descriptor,
                          var_type = "continuous_design"
                          )
        self.initial_point = initial_point
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.scale_type = scale_type
        self.scale_value = scale_value
class DiscreteDesignVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name)


class Simulation:
    def __init__(self,yaml_in,fname_out):
        self.hasHessiansImplemented  = False
        self.hasGradientsImplemented = False
        self.eval_id = {}
        self.num_parameters = {}
        self.num_functions  = {}
        self.parameters = {}
        self.asv = {}
        self.dvv = {}
        self.ac = {}
        self.return_values = {}
        self.readYamlString(yaml_in)

    def prerun(self):
        pass
        
    def run(self):
        self.return_values = {}
        for key in self.asv.keys():
            if (self.asv[key] & 1):
                self.calculateFunctions()
            if (self.asv[key] & 2):
                self.calculateGradients()
            if (self.asv[key] & 4):
                self.calculateHessians()

   def postrun(self):
       pass
    
    def readYamlString(self, yaml_string):
        dakota_in = yaml.load(yaml_string)
        self.eval_id       = dakota_in['eval_id']
        self.num_functions = dakota_in['num_functions']
        self.parameters    = dakota_in['variables']
        self.asv           = dakota_in['asv']
        self.dvv           = dakota_in['dvv']
        self.ac            = dakota_in['ac']

    def writeFile(self, fname_out):

        outfile = open(fname_out, 'w')
    
        # write functions
        keys = list(self.asv.keys())
        for func_ind in range(0, self.num_functions):
            if (self.asv[keys[func_ind]] & 1):
                functions = self.return_values['fns']    
                outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
        
        # write gradients
        for func_ind in range(0, self.num_functions):
            if (self.asv[keys[func_ind]] & 2):
                grad = self.return_values['fnGradients'][func_ind]
                outfile.write('[ ')
                for deriv in grad: 
                    outfile.write(str(deriv) + ' ')
                outfile.write(']\n')
        
        # write Hessians
        for func_ind in range(0, self.num_functions):
            if (self.asv[keys[func_ind]] & 4):
                hessian = self.return_values['fnHessians'][func_ind]
                outfile.write('[[ ')
                for hessrow in hessian:
                    for hesscol in hessrow:
                        outfile.write(str(hesscol) + ' ')
                    outfile.write('\n')
                outfile.write(']]')            

        outfile.close()


    def calculateFunctions(self):
        x0 = self.parameters['x1']
        x1 = self.parameters['x2']

        f0 = x1-x0*x0
        f1 = 1-x0
        
        f = [100*f0*f0+f1*f1]
        self.return_values['fns'] = f
        
    def calculateGradients(self):
        pass
        
    def calculateHessians(self):
        pass

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
