from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    exp,
    Expression,
    log,
    Set,
    Param,
    value,
    Var,
    units as pyunits,
    Reference,
    SolverFactory
)
from pyomo.common.config import ConfigBlock, ConfigDict, ConfigValue, In, Bool
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialFlowBasis,
    ElectrolytePropertySet,
)
from idaes.core.base.components import Component, __all_components__
from idaes.core.base.phases import (
    Phase,
    AqueousPhase,
    LiquidPhase,
    VaporPhase,
    __all_phases__,
)
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_activated_constraints,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.initialization.initializer_base import InitializerBase

from idaes.models.properties.modular_properties.base.generic_reaction import (
    equil_rxn_config,
)
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_phase_method,
    GenericPropertyPackageError,
    StateIndex,
    identify_VL_component_list,
    # estimate_Tbub,
    # estimate_Tdew,
    estimate_Pbub,
    estimate_Pdew,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType


def _init_Tbub(blk, T_units):
    for pp in blk.params._pe_pairs:
        l_phase, _, raoult_comps, henry_comps, _, _ = identify_VL_component_list(
            blk, pp
        )

        if raoult_comps == []:
            continue

        Tbub0 = _custom_estimate_Tbub(blk, T_units, raoult_comps, henry_comps, l_phase)

        blk.temperature_bubble[pp].set_value(Tbub0)

        for j in raoult_comps:
            blk._mole_frac_tbub[pp, j].value = value(
                blk.mole_frac_comp[j]
                * get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), Tbub0 * T_units
                )
                / blk.pressure
            )
            if blk.is_property_constructed("log_mole_frac_tbub"):
                blk.log_mole_frac_tbub[pp, j].value = value(
                    log(blk._mole_frac_tbub[pp, j])
                )

        for j in henry_comps:
            blk._mole_frac_tbub[pp, j].value = value(
                blk.mole_frac_comp[j]
                * blk.params.get_component(j)
                .config.henry_component[l_phase]["method"]
                .return_expression(blk, l_phase, j, Tbub0 * T_units)
                / blk.pressure
            )
            if blk.is_property_constructed("log_mole_frac_tbub"):
                blk.log_mole_frac_tbub[pp, j].value = value(
                    log(blk._mole_frac_tbub[pp, j])
                )


def _init_Tdew(blk, T_units):
    for pp in blk.params._pe_pairs:
        l_phase, _, raoult_comps, henry_comps, _, _ = identify_VL_component_list(
            blk, pp
        )

        if raoult_comps == []:
            continue

        Tdew0 = _custom_estimate_Tdew(blk, T_units, raoult_comps, henry_comps, l_phase)

        blk.temperature_dew[pp].set_value(Tdew0)

        for j in raoult_comps:
            blk._mole_frac_tdew[pp, j].value = value(
                blk.mole_frac_comp[j]
                * blk.pressure
                / get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), Tdew0 * T_units
                )
            )
            if blk.is_property_constructed("log_mole_frac_tdew"):
                blk.log_mole_frac_tdew[pp, j].value = value(
                    log(blk._mole_frac_tdew[pp, j])
                )
        for j in henry_comps:
            blk._mole_frac_tdew[pp, j].value = value(
                blk.mole_frac_comp[j]
                * blk.pressure
                / blk.params.get_component(j)
                .config.henry_component[l_phase]["method"]
                .return_expression(blk, l_phase, j, Tdew0 * T_units)
            )
            if blk.is_property_constructed("log_mole_frac_tdew"):
                blk.log_mole_frac_tdew[pp, j].value = value(
                    log(blk._mole_frac_tdew[pp, j])
                )


def _init_Pbub(blk):
    for pp in blk.params._pe_pairs:
        l_phase, _, raoult_comps, henry_comps, _, _ = identify_VL_component_list(
            blk, pp
        )

        if raoult_comps == []:
            continue

        blk.pressure_bubble[pp].set_value(
            estimate_Pbub(blk, raoult_comps, henry_comps, l_phase)
        )

        for j in raoult_comps:
            blk._mole_frac_pbub[pp, j].value = value(
                blk.mole_frac_comp[j]
                * blk.pressure_sat_comp[j]
                / blk.pressure_bubble[pp]
            )
            if blk.is_property_constructed("log_mole_frac_pbub"):
                blk.log_mole_frac_pbub[pp, j].value = value(
                    log(blk._mole_frac_pbub[pp, j])
                )
        for j in henry_comps:
            blk._mole_frac_pbub[pp, j].value = value(
                blk.mole_frac_comp[j] * blk.henry[l_phase, j] / blk.pressure_bubble[pp]
            )
            if blk.is_property_constructed("log_mole_frac_pbub"):
                blk.log_mole_frac_pbub[pp, j].value = value(
                    log(blk._mole_frac_pbub[pp, j])
                )


def _init_Pdew(blk):
    for pp in blk.params._pe_pairs:
        l_phase, _, raoult_comps, henry_comps, _, _ = identify_VL_component_list(
            blk, pp
        )

        if raoult_comps == []:
            continue

        blk.pressure_dew[pp].set_value(
            estimate_Pdew(blk, raoult_comps, henry_comps, l_phase)
        )

        for j in raoult_comps:
            blk._mole_frac_pdew[pp, j].value = value(
                blk.mole_frac_comp[j] * blk.pressure_dew[pp] / blk.pressure_sat_comp[j]
            )
            if blk.is_property_constructed("log_mole_frac_pdew"):
                blk.log_mole_frac_pdew[pp, j].value = value(
                    log(blk._mole_frac_pdew[pp, j])
                )
        for j in henry_comps:
            blk._mole_frac_pdew[pp, j].value = value(
                blk.mole_frac_comp[j] * blk.pressure_dew[pp] / blk.henry[l_phase, j]
            )
            if blk.is_property_constructed("log_mole_frac_pdew"):
                blk.log_mole_frac_pdew[pp, j].value = value(
                    log(blk._mole_frac_pdew[pp, j])
                )



from enum import Enum

from pyomo.environ import units as pyunits, value

from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
)
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


TOL = 1e-1
MAX_ITER = 30


def _custom_estimate_Tbub(blk, T_units, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate bubble point temperature

    Args:
        blk: StateBlock to use
        T_units: units of temperature
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated bubble point temperature as a float.

    """
    # Use lowest component temperature_crit as starting point
    # Starting high and moving down generally works better,
    # as it under-predicts next step due to exponential form of
    # Psat.
    # Subtract 1 to avoid potential singularities at Tcrit
    Tbub_initial = (
        min(blk.params.get_component(j).temperature_crit.value for j in raoult_comps)
        - 1
    )

    m = Block()
    blk.add_component("tbub_initialise", m)
    m.Tbub = Var(initialize=Tbub_initial)
    m.f = Expression(expr=(
        sum(
            get_method(blk, "pressure_sat_comp", j)(
                blk, blk.params.get_component(j), m.Tbub * T_units
            )
            * blk.mole_frac_comp[j]
            for j in raoult_comps
        )
        + sum(
            blk.mole_frac_comp[j]
            * blk.params.get_component(j)
            .config.henry_component[liquid_phase]["method"]
            .return_expression(blk, liquid_phase, j, m.Tbub * T_units)
            for j in henry_comps
        )
        - blk.pressure
    ))
    m.df = Expression(expr=(
        sum(
            get_method(blk, "pressure_sat_comp", j)(
                blk, blk.params.get_component(j), m.Tbub * T_units, dT=True
            )
            * blk.mole_frac_comp[j]
            for j in raoult_comps
        )
        + sum(
            blk.mole_frac_comp[j]
            * blk.params.get_component(j)
            .config.henry_component[liquid_phase]["method"]
            .dT_expression(blk, liquid_phase, j, m.Tbub * T_units)
            for j in henry_comps
        )
    ))
    m.tbub_constraint = Constraint(expr=(
        m.Tbub == m.Tbub - m.f / m.df
    ))
    solver = SolverFactory("ipopt")
    solver.options["tol"] = 1e-8
    solver.solve(m)

    result = m.Tbub.value
    blk.del_component(m)  # cleanup
    return result


def _custom_estimate_Tdew(blk, T_units, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate dew point temperature

    Args:
        blk: StateBlock to use
        T_units: units of temperature
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated dew point temperature as a float.

    """
    # Use lowest component critical temperature
    # as starting point
    # Subtract 1 to avoid potential singularities at Tcrit
    Tdew_initial = (
        min(blk.params.get_component(j).temperature_crit.value for j in raoult_comps)
        - 1
    )

    m = Block()
    blk.add_component("tdew_initialise", m)
    m.Tdew = Var(initialize=Tdew_initial)
    m.f = Expression(expr=(
        blk.pressure
        * (
            sum(
                blk.mole_frac_comp[j]
                / get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), m.Tdew * T_units
                )
                for j in raoult_comps
            )
            + sum(
                blk.mole_frac_comp[j]
                / blk.params.get_component(j)
                .config.henry_component[liquid_phase]["method"]
                .return_expression(blk, liquid_phase, j, m.Tdew * T_units)
                for j in henry_comps
            )
        )
        - 1
    ))
    m.df = Expression(expr=-(
        blk.pressure
        * (
            sum(
                blk.mole_frac_comp[j]
                / get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), m.Tdew * T_units
                )
                ** 2
                * get_method(blk, "pressure_sat_comp", j)(
                    blk,
                    blk.params.get_component(j),
                    m.Tdew * T_units,
                    dT=True,
                )
                for j in raoult_comps
            )
            + sum(
                blk.mole_frac_comp[j]
                / blk.params.get_component(j)
                .config.henry_component[liquid_phase]["method"]
                .return_expression(blk, liquid_phase, j, m.Tdew * T_units)
                ** 2
                * blk.params.get_component(j)
                .config.henry_component[liquid_phase]["method"]
                .dT_expression(blk, liquid_phase, j, m,Tdew * T_units)
                for j in henry_comps
            )
        )
    ))
    m.Tdew_constraint = Constraint(expr=(
        m.Tdew == m.Tdew - m.f / m.df
    ))
    solver = SolverFactory("ipopt")
    solver.options["tol"] = TOL
    solver.solve(m)

    result = m.Tdew.value
    blk.del_component(m)  # cleanup
    return result
