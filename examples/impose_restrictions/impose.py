# -*- coding: utf-8 -*

class Constraint:
    def __init__(self,var_name,condition,var_value):            
        self.var_name = var_name
        self.condition = condition
        self.value = value
        
        supported_conditions = ['<','>','=','!=', '<=', '>=']
        
        if self.condition not in supported_conditions:
            err_msg = "{} is not a supported equality inequality condition"
            raise ValueError(err_msg)
        
    def __validate(self,condition):
        pass

class ConstraintList:
    def __init__(self):
        self.constraints = []
        
    def add_constraint(self,
                       var_name, 
                       condition, 
                       var_value):
        
        """
        Arguments:
        
        var_name (string) - variable name
        var_name (string) - condition
        var_name (string) - variable value
        """
        self.constraints.append(Constraint(var_name,
                                           condition,
                                           var_value))
    
    def string(self):
        str_out = ""
        for c in self.constraints:
            str_out += "{} {} {}".format(c.var_name,
                                         c.condition,
                                         c.value)
        return str_out
                                         