from pyomo.environ import (
    Expression,
    SolverFactory,
    check_optimal_termination,
    units,
    Block,
    Constraint,
    value
)
from idaes.models.properties.general_helmholtz.helmholtz_state import HelmholtzStateBlockData, _StateBlock
from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzParameterBlockData
from idaes.models.properties.general_helmholtz import HelmholtzThermoExpressions
from idaes.core import declare_process_block_class
from idaes.core.util.initialization import solve_indexed_blocks, revert_state_vars
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from ahuora_property_packages.base.state_block_constraints import StateBlockConstraints
from ahuora_property_packages.utils.fix_state_vars import fix_state_vars


class _ExtendedStateBlock(_StateBlock):

    def initialize(blk, *args, **kwargs):
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        set_vapor_frac_guesses(blk)

        flag_dict = fix_state_vars(blk, kwargs.get("state_args", None))

        dof = degrees_of_freedom(blk)
        if dof != 0:
            raise InitializationError(
                f"{blk.name} Unexpected degrees of freedom during "
                f"initialization at property initialization step: {dof}."
            )
        
        res = None
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
        
        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        if kwargs.get("hold_state") is True:
            return flag_dict
        else:
            blk.release_state(flag_dict)
        
        init_log.info(
            "Property package initialization: {}.".format(idaeslog.condition(res))
        )
    
    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        revert_state_vars(blk, flags)



def set_vapor_frac_guesses(blk: Block) -> None:
    """
    Vapor fraction specified
    valid combinations: {p, x}, {T, x}

    We find a guess for h and p based on the specified vapor fraction and either p or T.

    This is because the initial guess has to be within the vapor fraction region
    otherwise the solver has a hard time converging.
    """
    helmholtz_blk = blk[blk.index_set().first()].config.parameters
    
    for sb in blk.values():
        if not hasattr(sb.constraints, "vapor_frac"):
            continue
        
        # Get the value the vapor fraction is set to.
        # value(sb.constraints.vapor_frac) returns the current value of the left hand
        # side of the constraint, which (since the model is not solved) is not actually
        # the correct value of the vapor fraction.
        # instead, the lower and upper bounds of the constraint are equal and set to the 
        # desired value of the vapor fraction.
        if sb.constraints.vapor_frac.lower != sb.constraints.vapor_frac.upper:
            raise ValueError("We haven't implemented supporting bounds on the vapor fraction constraint yet. See PropertyPackages, helmholtz_extended.py")
        x = sb.constraints.vapor_frac.lower
        del sb.constraints.vapor_frac # new smooth constraint will be added later instead
        

        # add custom vapor fraction constraint: h = h_sat_liq + x(h_sat_vap - h_sat_liq)
        # this is to make the vapor fraction continuous, rather than cutting off at 0 and 1. Helps with solving reliability as outside
        # the vapor fraction region, this will still be smooth and will not have a zero gradient.
        sb.constraints.add_component(
            "custom_vapor_frac",
            Constraint(expr=sb.enth_mol == sb.enth_mol_sat_phase["Liq"] + x * (sb.enth_mol_sat_phase["Vap"] - sb.enth_mol_sat_phase["Liq"]))
        )

@declare_process_block_class("HelmholtzExtendedStateBlock", block_class=_ExtendedStateBlock)
class HelmholtzExtendedStateBlockData(HelmholtzStateBlockData, StateBlockConstraints):

    def build(blk, *args):
        HelmholtzStateBlockData.build(blk, *args)
        StateBlockConstraints.build(blk, *args)

    def add_extra_expressions(blk):
        super().add_extra_expressions()

        # Generic property packages support temperature_bubble and temperature_dew,
        # Helmholtz is pure so these are the same, and so it only has
        # temperature_sat. https://idaes-pse.readthedocs.io/en/stable/reference_guides/model_libraries/generic/property_models/helmholtz.html#expressions
        # We add them here so that the property package interface is consistent.
        blk.add_component("temperature_bubble", Expression(expr=blk.temperature_sat))
        blk.add_component("temperature_dew", Expression(expr=blk.temperature_sat))

        # For numerical stability when constraining temperature, we actually smooth out the temperature equilibrium curve. 
        # Otherwise, the solver has difficulty solving through the vapor fraction region as there is zero gradient.
    
        # rename temperature to old_temperature
        old_temperature = blk.temperature
        blk.del_component("temperature")
        blk.add_component("old_temperature", old_temperature)
        # create smooth temperature expression
        blk.add_component("temperature", Expression(
            expr=blk.old_temperature + (blk.enth_mol / (1 * units.J/units.mol))*0.000001 * units.K
        ))


@declare_process_block_class("HelmholtzExtendedParameterBlock")
class HelmholtzExtendedParameterBlockData(HelmholtzParameterBlockData):
    def build(self):
        super().build()
        self._state_block_class = HelmholtzExtendedStateBlock # noqa: F821
