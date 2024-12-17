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
from idaes.core.util.model_statistics import degrees_of_freedom
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
        print("------------------------------------")
        print("Starting initialization of ", self, "dof", degrees_of_freedom(self)) 
        print("------------------------------------")

        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        if (degrees_of_freedom(self) == 0):
            # TODO: this per block, rather than degrees of freedom for the whole model
            print("trying initial solve")
            for i,b in self.items():
                print(degrees_of_freedom(b))
                s = SolverFactory('ipopt')
                res = s.solve(b, tee=False)
                print(res.solver.termination_condition)
            #return {} # Trying just to still run normal initialization anyways

        print("Initial State")
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")


        
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



        print("State at deactivation")
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("Dof before:", degrees_of_freedom(self))
        print("-----")

        res = super().initialize(*args, **kwargs)
        print("dof after initialise:" , degrees_of_freedom(self))
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("-----")
        super().release_state(res) # release the state after initialization, so we can just fix the variables we need to fix.
        print("Dof after release:",degrees_of_freedom(self))

        print("State after initialization")
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("-----")



        flags = {}
        [print(v) for v in self.items()]
        # reactivate the variables that were deactivated
        for i, b in self.items():
            for name, value, in deactivated_vars[i].items():
                print("reactivating",name)
                getattr(b,name).fix(value)

        # if zero degrees of freedom, resovle the model with these constraints
        # Todo: resolve the model with whatever constraints are active, and guess the rest if degrees of freedom > 0
        # for i,b in self.items():
        #     if degrees_of_freedom(b) == 0:
        #         s = SolverFactory('ipopt')
        #         res = s.solve(b, tee=True)
        #         print("init:11o",res.solver.termination_condition)

        for i, b in self.items():
            print("activating constraints")
            b.constraints.activate()
            flags[i] = {}
            if hold_state:
                # Fix the required state variables for zero degrees of freedom, and return a dictionary of the flags.
                if not hasattr(b.constraints,"flow_mass") and not b.flow_mol.is_fixed():
                    # We need to fix the flow_mol variable
                    flags[i]["flow_mol"] = True
                    b.flow_mol.fix()
                avaliable_constraints = ["enth_mass","entr_mass","temperature","entr_mol","smooth_temperature","vapor_frac"]
                if not b.enth_mol.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for prop_name in avaliable_constraints:
                        # check if it's a constraint that 
                        if hasattr(b.constraints,prop_name):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(prop_name)
                            found_constraint = True
                            break
                        if hasattr(b,prop_name):
                            # check if it's a variable that needs to be unfixed
                            if type(getattr(b,prop_name)) in (ScalarVar,):
                                var = getattr(b,prop_name)
                                if var.is_fixed():
                                    # we need to move this from the list of constraints,
                                    # it's used to define the state
                                    avaliable_constraints.remove(prop_name)
                                    found_constraint = True
                                    break


                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["enth_mol"] = True
                        print("fixing enthalpy")
                        b.enth_mol.fix()
                if not b.pressure.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for prop_name in avaliable_constraints:
                        if hasattr(b.constraints,prop_name):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(prop_name)
                            found_constraint = True
                            break
                        if hasattr(b,prop_name):
                            # check if it's a variable that needs to be unfixed
                            if type(getattr(b,prop_name)) in (ScalarVar,):
                                var = getattr(b,prop_name)
                                if var.is_fixed():
                                    # we need to move this from the list of constraints,
                                    # it's used to define the state
                                    avaliable_constraints.remove(prop_name)
                                    found_constraint = True
                                    break
                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["pressure"] = True
                        b.pressure.fix() 
        # Solve again with new constraints
        # print("State at end initialization,hold state = ", hold_state,"dof", degrees_of_freedom(self))
        # for i,b in self.items():
        #     if degrees_of_freedom(b) == 0:
        #         print("solving again")
        #         s = SolverFactory('ipopt')
        #         res = s.solve(b, tee=True)
        #         print("init:end",res.solver.termination_condition)
        # for i,b in self.items():
        #     for var in b.component_data_objects(Var):
        #         if var.is_fixed():
        #             print(f"{var.name} = {var.value}")
        print("-----")
        print("Flags",flags)

        print("------------------------------------")
        print("init state end", self)
        print(degrees_of_freedom(self))
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("------------------------------------")
        return flags
    
    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        print("------------------------------------")
        print("release state,", self, flags)
        print(degrees_of_freedom(self))
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("------------------------------------")

        
        print([b.name for i,b in self.items()])
        for i, v in self.items():
            if i not in flags:
                continue
            for key in flags[i]:
                getattr(v,key).unfix()
        print("release state end")
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        # REsolve block
        # for i,b in self.items():
        #     if degrees_of_freedom(b) == 0:
        #         print("release_state solving again")
        #         s = SolverFactory('ipopt')
        #         res = s.solve(b, tee=True)
        #         print("release:end", b.name ,res.solver.termination_condition)
        print("------------------------------------")
        print("release state end", self)
        print(degrees_of_freedom(self))
        for i,b in self.items():
            for var in b.component_data_objects(Var):
                if var.is_fixed():
                    print(f"{var.name} = {var.value}")
        print("------------------------------------")

    

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
