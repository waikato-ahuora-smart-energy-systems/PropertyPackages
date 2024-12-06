# extends the IDAES helmholtz property package to include additional properties and methods.
from idaes.models.properties.general_helmholtz.helmholtz_state import HelmholtzStateBlockData, _StateBlock
from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzParameterBlockData
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock, _GenericStateBlock, GenericParameterData, GenericStateBlockData
from idaes.core import declare_process_block_class
from property_packages.utils.add_extra_expressions import add_extra_expressions
from pyomo.environ import Constraint, Block
from pyomo.core.base.expression import Expression, ScalarExpression, _GeneralExpressionData, ExpressionData
from pyomo.core.base.var import IndexedVar, ScalarVar, Var, _GeneralVarData,VarData
import idaes.logger as idaeslog
from pyomo.environ import SolverFactory

# NOTE:
# THis only works for FTPx formulation right now.




class _ExtendedGenericStateBlock(_GenericStateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # Missing argument

    def initialize(self, *args, **kwargs):
        hold_state = kwargs.pop("hold_state", False)
        deactivated_vars = {}
        for i, b in self.items():
            b.constraints.deactivate()
            # deactivate anything that is not a state var, because they'll ruin the initialization
            # iterate over all vars in block 
            state_vars = b.define_state_vars()
            deactivated_vars[i] = {}
            for var in b.component_data_objects(Var, descend_into=False):
                if var.local_name not in state_vars and type(var) in (ScalarVar,): # TODO: Support this for indexed vars e.g mole_frac
                    if var.is_fixed():
                        deactivated_vars[i][var.local_name] = var.value
                        var.unfix()
                    
                
        kwargs["hold_state"] = True # this is done to avoid release_state being called by super().initialize
        res = super().initialize(*args, **kwargs)
        super().release_state(res) # release the state after initialization, so we can just fix the variables we need to fix.
        flags = {}
        print(self.items())
        # reactivate the variables that were deactivated
        for i, b in self.items():
            for name, value, in deactivated_vars[i].items():
                getattr(b,name).fix(value)

        for i, b in self.items():
            b.constraints.activate()
            flags[i] = {}
            if hold_state:
                # Fix the required state variables for zero degrees of freedom, and return a dictionary of the flags.
                if not hasattr(b.constraints,"flow_mass") and not b.flow_mol.is_fixed():
                    # We need to fix the flow_mol variable
                    flags[i]["flow_mol"] = True
                    b.flow_mol.fix()
                avaliable_constraints = ["enth_mass","temperature","entr_mass","entr_mol","smooth_temperature","vapor_frac"]
                if not b.enth_mol.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for constraint in avaliable_constraints:
                        if hasattr(b.constraints,constraint):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(constraint)
                            found_constraint = True
                            break
                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["enth_mol"] = True
                        b.enth_mol.fix()
                if not b.pressure.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for constraint in avaliable_constraints:
                        if hasattr(b.constraints,constraint):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(constraint)
                            found_constraint = True
                            break
                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["pressure"] = True
                        b.pressure.fix() 
        return flags
    
    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        for i, v in self.items():
            for key in flags[i]:
                getattr(v,key).unfix()
    

@declare_process_block_class("GenericExtendedStateBlock", block_class=_ExtendedGenericStateBlock)
class GenericExtendedStateBlockData(GenericStateBlockData):

    def build(self, *args):
        super().build(*args)
        # Add expressions for smooth_temperature, enthalpy in terms of mass, etc.
        add_extra_expressions(self)
        # Add a block for constraints, so we can disable or enable them in bulk
        self.constraints = Block()
    
    def constrain(self,name:str,value:float):
        # Value must be a float. TODO: Handle unit conversion.
        var = getattr(self,name)
        if type(var) == ScalarExpression:
            self.constraints.add_component(name, Constraint(expr=var == value))
        elif type(var) in (ScalarVar, _GeneralVarData, VarData):
            var.fix(value)
        elif type(var) in ( _GeneralExpressionData, ExpressionData) :
            # allowed, but we don't need to fix it (eg. mole_frac_comp in helmholtz)
            print(f"Variable {self} {name} is an Expression: {type(var)}")
            pass
        else:
            raise Exception(f"Variable {self} {name} is not a Var or Expression: {type(var)}")

@declare_process_block_class("GenericExtendedParameterBlock")
class GenericExtendedParameterData(GenericParameterData):
    def build(self):
        super().build()
        # set that we should use the extended state block
        self._state_block_class = GenericExtendedStateBlock # type: ignore because it'll get created with declare_process_block_class
