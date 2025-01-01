from pyomo.environ import Constraint, Block, SolverFactory
from pyomo.core.base.var import IndexedVar, ScalarVar, Var, _GeneralVarData, VarData
from pyomo.core.base.expression import Expression, ScalarExpression, _GeneralExpressionData, ExpressionData

from idaes.models.properties.general_helmholtz.helmholtz_state import HelmholtzStateBlockData, _StateBlock
from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzParameterBlockData
from idaes.core import declare_process_block_class
from idaes.core.util.initialization import solve_indexed_blocks, revert_state_vars
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from property_packages.utils.add_extra_expressions import add_extra_expressions


def fix_state_vars(blk, state_args=None):
    """
    Method for fixing state variables within StateBlocks. Method takes an
    optional argument of values to use when fixing variables.

    This method is analogous to the fix_state_vars method in
    idaes.core.util.initialization, but it allows us to handle
    the case where we have constraints that define the state.

    Args:
        blk : An IDAES StateBlock object in which to fix the state variables
        state_args : a dict containing values to use when fixing state
                variables. Keys must match with names used in the
                define_state_vars method, and indices of any variables must
                agree.

    Returns:
        A dict keyed by block index, state variable name (as defined by
        define_state_variables) and variable index indicating the fixed status
        of each variable before the fix_state_vars method was applied.
    """
    if state_args is None:
        state_args = {}

    flags = {}
    for k, b in blk.items():
        available_constraints = [
            "temperature",
            "smooth_temperature",
            "enth_mass",
            "entr_mol",
            "entr_mass",
            "total_energy_flow",
            "custom_vapor_frac",
        ]
        for n, v in b.define_state_vars().items():
            fix_var = True
            if n == "flow_mol":
                for prop_name in ["flow_mass", "flow_vol"]:
                    if hasattr(b.constraints, prop_name):
                        fix_var = False
                        break
            elif n in ["enth_mol", "pressure"]:
                if not v.is_fixed():
                    # check if any of the constraints exist
                    for prop_name in available_constraints:
                        if hasattr(b.constraints, prop_name):
                            # don't fix this variable - it is defined by a constraint
                            available_constraints.remove(prop_name)
                            fix_var = False
                            break
            for i in v:
                flags[k, n, i] = v[i].is_fixed()
                if fix_var and not v[i].is_fixed():
                    # If not fixed, fix at either guess provided or current value
                    if n in state_args:
                        # Try to get initial guess from state_args
                        try:
                            if i is None:
                                val = state_args[n]
                            else:
                                val = state_args[n][i]
                        except KeyError:
                            raise ConfigurationError(
                                "Indexes in state_args did not agree with "
                                "those of state variable {}. Please ensure "
                                "that indexes for initial guesses are correct.".format(
                                    n
                                )
                            )
                        v[i].fix(val)
                    else:
                        # No guess, try to use current value
                        if v[i].value is not None:
                            v[i].fix()
                        else:
                            # No initial value - raise Exception before this
                            # gets to a solver.
                            raise ConfigurationError(
                                "State variable {} does not have a value "
                                "assigned. This usually occurs when a Var "
                                "is not assigned an initial value when it is "
                                "created. Please ensure all variables have "
                                "valid values before fixing them.".format(v.name)
                            )
    return flags


class _ExtendedStateBlock(_StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) 

    def initialize(blk, *args, **kwargs):
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        flag_dict = fix_state_vars(blk, kwargs.get("state_args", None))

        dof = degrees_of_freedom(blk)
        if dof != 0:
            raise InitializationError(
                f"{blk.name} Unexpected degrees of freedom during "
                f"initialization at property initialization step: {dof}."
            )
        opt = SolverFactory('ipopt')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            try:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            except ValueError as e:
                if str(e).startswith("No variables appear"):
                    # https://github.com/Pyomo/pyomo/pull/3445
                    pass
                else:
                    raise e

        if kwargs.get("hold_state") is True:
            return flag_dict
        else:
            blk.release_state(flag_dict)
    
    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        revert_state_vars(blk, flags)


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
