# extends the IDAES helmholtz property package to include additional properties and methods.
from idaes.models.properties.general_helmholtz.helmholtz_state import HelmholtzStateBlockData, _StateBlock
from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzParameterBlockData
from idaes.core import declare_process_block_class
from property_packages.utils.add_extra_expressions import add_extra_expressions
from pyomo.environ import Constraint, Block
from pyomo.core.base.expression import Expression, ScalarExpression, _GeneralExpressionData, ExpressionData
from pyomo.core.base.var import IndexedVar, ScalarVar, Var, _GeneralVarData,VarData
import idaes.logger as idaeslog


class _ExtendedStateBlock(_StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) 

    def initialize(self, *args, **kwargs):
        hold_state = kwargs.pop("hold_state", False)
        for i, v in self.items():
            v.constraints.deactivate()
        res = super().initialize(*args, **kwargs)
        flags = {}
        for i, v in self.items():
            v.constraints.activate()
            flags[i] = {}
            if hold_state:
                # Fix the required state variables for zero degrees of freedom, and return a dictionary of the flags.
                if not hasattr(v.constraints,"flow_mass") and not v.flow_mol.is_fixed():
                    # We need to fix the flow_mol variable
                    flags[i]["flow_mol"] = True
                    v.flow_mol.fix()
                avaliable_constraints = ["enth_mass","temperature","total_energy_flow","entr_mass","entr_mol","smooth_temperature","custom_vapor_frac"]
                if not v.enth_mol.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for constraint in avaliable_constraints:
                        if hasattr(v.constraints,constraint):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(constraint)
                            found_constraint = True
                            break
                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["enth_mol"] = True
                        v.enth_mol.fix()
                if not v.pressure.is_fixed():
                    # check if any of the constraints exist
                    found_constraint = False
                    for constraint in avaliable_constraints:
                        if hasattr(v.constraints,constraint):
                            # we don't need to fix the variable
                            # but we need to remove this from the list of constraints (it can't be used to define pressure)
                            avaliable_constraints.remove(constraint)
                            found_constraint = True
                            break
                    if not found_constraint:
                        # we need to fix the variable
                        flags[i]["pressure"] = True
                        v.pressure.fix() 
        return flags
    
    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        for i, v in self.items():
            for key in flags[i]:
                getattr(v,key).unfix()
    

@declare_process_block_class("HelmholtzExtendedStateBlock", block_class=_ExtendedStateBlock)
class HelmholtzExtendedStateBlockData(HelmholtzStateBlockData):

    def build(self, *args):
        super().build(*args)
        # Add expressions for smooth_temperature, enthalpy in terms of mass, etc.
        add_extra_expressions(self)
        # Add a block for constraints, so we can disable or enable them in bulk
        self.constraints = Block()


    def constrain(self, name: str, value: float) -> Constraint | Var | None:
        """constrain a component by name to a value"""
        # TODO: handle unit conversion
        var = getattr(self, name)
        return self.constrain_component(var, value)


    def constrain_component(self, component: Var | Expression, value: float) -> Constraint | Var | None:
        """
        Constrain a component to a value
        """
        if type(component) == ScalarExpression:
            c = Constraint(expr=component == value)
            self.constraints.add_component(component.local_name, c)
            return c
        elif type(component) in (ScalarVar, _GeneralVarData, VarData, IndexedVar):
            component.fix(value)
            return component
        elif type(component) in (_GeneralExpressionData, ExpressionData):
            # allowed, but we don't need to fix it (eg. mole_frac_comp in helmholtz)
            return None
        else:
            raise Exception(
                f"Component {component} is not a Var or Expression: {type(component)}"
            )

@declare_process_block_class("HelmholtzExtendedParameterBlock")
class HelmholtzExtendedParameterBlockData(HelmholtzParameterBlockData):
    def build(self):
        super().build()
        self._state_block_class = HelmholtzExtendedStateBlock # noqa: F821
