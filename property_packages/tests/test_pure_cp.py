# Build and solve a state block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock, Component
from idaes.core.solvers import get_solver

from CoolProp import CoolProp

from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import report_statistics, degrees_of_freedom
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType as PT
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
import idaes.core.util.scaling as iscale
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant

from math import floor
from typing import Any, Dict, List
from compounds.Compound import Compound
from pyomo.environ import units as pyunits
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (LogBubbleDew)
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType as PT
from idaes.models.properties.modular_properties.phase_equil import (SmoothVLE)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.pure import RPP4, RPP3, Perrys
from property_packages.modular.builder.data.chem_sep import ChemSep
from pyomo.common.fileutils import this_file_dir
from property_packages.types import States
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
import csv
from numpy import arange
from pyomo.environ import (
    Block,
    ConcreteModel,
    Param,
    SolverStatus,
    TerminationCondition,
    units as pyunits,
    value,
    Var,
)
import pytest

solver = get_solver(solver="ipopt_v2")

from idaes.models.properties.modular_properties.coolprop.coolprop_wrapper import (
    CoolPropWrapper,
    CoolPropExpressionError,
    CoolPropPropertyError,
)
from CoolProp.CoolProp import PropsSI
import idaes.logger as idaeslog
SOUT = idaeslog.INFO

def get_m(name="benzene"):

    # Clear cached components to ensure clean slate
    CoolPropWrapper.flush_cached_components()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = build_package("coolprop", [name], ["Liq"])
    m.fs.state = m.fs.props.build_state_block([0], defined_state=True)
    m.fs.state[0].flow_mol.fix(1)
    m.fs.state[0].mole_frac_comp[name].fix(1)
    return m

# def test_cubic_liquid():
#     m = get_m()

#     for P in range(1, 11):
#         Td = CoolProp.PropsSI("T", "P", P * 1e5, "Q", 0.5, "PR::benzene")
#         Tmin = CoolProp.PropsSI("TMIN", "benzene")

#         m.fs.state[0].pressure.fix(P * 1e5)
#         m.fs.state[0].temperature.fix(Tmin)

#         m.fs.state.initialize()

#         for T in arange(Tmin, Td, 10):
#             m.fs.state[0].temperature.fix(T)

#             results = solver.solve(m.fs)

#             assert (
#                 results.solver.termination_condition == TerminationCondition.optimal
#             )
#             assert results.solver.status == SolverStatus.ok

#             # Check results
#             assert pytest.approx(
#                 CoolProp.PropsSI("Z", "T", T, "P", P * 1e5, "PR::benzene"), rel=1e-8
#             ) == value(m.fs.state[0].compress_fact_phase["Liq"])

#             assert pytest.approx(
#                 CoolProp.PropsSI("DMOLAR", "T", T, "P", P * 1e5, "PR::benzene"),
#                 rel=1e-6,
#             ) == value(m.fs.state[0].dens_mol_phase["Liq"])

#             assert pytest.approx(
#                 CoolProp.PropsSI(
#                     "HMOLAR_RESIDUAL", "T", T, "P", P * 1e5, "PR::benzene"
#                 ),
#                 rel=1e-6,
#             ) == value(m.fs.state[0].enth_mol_phase["Liq"])

# def test_cubic_liquid_entr():
#     m = get_m()

#     Td = CoolProp.PropsSI("T", "P", 101325, "Q", 0.5, "PR::benzene")
#     Tmin = CoolProp.PropsSI("TMIN", "benzene")

#     for T in arange(Tmin, Td, 10):
#         m.fs.state[0].pressure.fix(101325)
#         m.fs.state[0].temperature.fix(T)

#         m.fs.state.initialize()

#         S0_CP = CoolProp.PropsSI("SMOLAR", "T", T, "P", 101325, "PR::benzene")
#         S0_I = value(m.fs.state[0].entr_mol_phase["Liq"])

#         for P in range(1, 11):
#             m.fs.state[0].pressure.fix(P * 1e5)

#             results = solver.solve(m.fs)

#             assert (
#                 results.solver.termination_condition == TerminationCondition.optimal
#             )
#             assert results.solver.status == SolverStatus.ok

#             assert pytest.approx(
#                 CoolProp.PropsSI("SMOLAR", "T", T, "P", P * 1e5, "PR::benzene") - S0_CP,
#                 rel=1e-4,
#             ) == value(m.fs.state[0].entr_mol_phase["Liq"] - S0_I)

# def test_cubic_liquid_enth():

#     m = get_m()

#     Td = CoolProp.PropsSI("T", "P", 101325, "Q", 0.5, "PR::benzene")
#     Tmin = CoolProp.PropsSI("TMIN", "benzene")

#     for T in arange(Tmin, Td, 10):
#         m.fs.state[0].pressure.fix(101325)
#         m.fs.state[0].temperature.fix(T)

#         m.fs.state.initialize()

#         S0_CP = CoolProp.PropsSI("HMOLAR", "T", T, "P", 101325, "PR::benzene")
#         S0_I = value(m.fs.state[0].enth_mol_phase["Liq"])

#         for P in range(1, 11):
#             m.fs.state[0].pressure.fix(P * 1e5)

#             results = solver.solve(m.fs)

#             assert (
#                 results.solver.termination_condition == TerminationCondition.optimal
#             )
#             assert results.solver.status == SolverStatus.ok

#             assert pytest.approx(
#                 CoolProp.PropsSI("HMOLAR", "T", T, "P", P * 1e5, "PR::benzene") - S0_CP,
#                 rel=1e-4,
#             ) == value(m.fs.state[0].enth_mol_phase["Liq"] - S0_I)


def test_compounds():

    for compound in ["benzene", "toluene", "water", "carbon dioxide"]:
        cool_prop_tester(compound)


def cool_prop_tester(name):
    """
    Method tests build_package("coolprop", [name], ["Liq"]) against CoolProp.PropsSI values

    Parameters:
        name (str): Name of cool-prop compound to be tested e.g. "benzene"

    Returns:
        none: Method asserts values against CoolProp.PropsSI
    """

    m = get_m(name)

    for P in range(1, 11):

        Td = CoolProp.PropsSI("T", "P", P * 1e5, "Q", 0.5, f"PR::{name}")
        Tmin = CoolProp.PropsSI("TMIN", name)

        m.fs.state[0].pressure.fix(P * 1e5)
        m.fs.state[0].temperature.fix(Tmin)

        m.fs.state.initialize()

        H0_CP = CoolProp.PropsSI("HMOLAR", "T", Tmin, "P", P*1e5, f"PR::{name}")
        H0_I = value(m.fs.state[0].enth_mol_phase["Liq"])

        S0_CP = CoolProp.PropsSI("SMOLAR", "T", Tmin, "P", P*1e5, f"PR::{name}")
        S0_I = value(m.fs.state[0].entr_mol_phase["Liq"])

        for T in arange(Tmin, Td, 10):
            m.fs.state[0].temperature.fix(T)

            results = solver.solve(m.fs)

            # Ensuring model solves optimally
            assert results.solver.termination_condition == TerminationCondition.optimal
            assert results.solver.status == SolverStatus.ok

            # Check results
            
            # Molar Density
            assert PropsSI("DMOLAR", "T", T, "P", P * 1e5, f"PR::{name}") == approx(value(m.fs.state[0].dens_mol_phase["Liq"])), f"Molar Density failed"
            # Molar Enthalpy
            assert PropsSI("HMOLAR", "T", T, "P", P * 1e5, f"PR::{name}") - H0_CP == approx(value(m.fs.state[0].enth_mol_phase["Liq"]) - H0_I)
            # Molar Entropy
            assert PropsSI("SMOLAR", "T", T, "P", P * 1e5, f"PR::{name}") - S0_CP == approx(value(m.fs.state[0].entr_mol_phase["Liq"]) - S0_I), f"{P, T, name}"

            # Internal Energy
            assert PropsSI("UMOLAR", "T", T, "P", P * 1e5, f"PR::{name}") == approx(value(m.fs.state[0].int_energy_mol_phase["Liq"])), f"Internal Energy failed"
            # Acentric Factor
            assert PropsSI("ACENTRIC", name) == approx(value(m.fs.props.acentric_factor[name]))

            # -- Ideal Helmholtz energy --
            # Specific Heat Capacity (cp)
            assert PropsSI("CP0MOLAR", "T", T, "P", P * 1e5, f"PR::{name}") == approx(value(m.fs.state[0].cp_mol_ig_comp))
            # Specific Heat Capacity (cv) <- doesn't exist
            # -- Virial EOS --
            # Molar gas constant
            assert PropsSI("GAS_CONSTANT", name) == approx(value(m.fs.props.gas_const))
            # Critical Pressurre
            assert PropsSI("P_CRITICAL", name) == approx(value(m.fs.props.pressure_crit[name]))
            # Phase Index
            assert PropsSI("PHASE_INDEX", name) == approx(value(m.fs.props.phase[name]))
            # Triple Point Pressure
            assert PropsSI("P_TRIPLE", name) == approx(value(m.fs.props.pressure_triple[name]))
            # Critical Temperature
            assert PropsSI("T_CRITICAL", name) == approx(value(m.fs.props.temperature_crit[name]))
            # Triple Point Temperature
            assert PropsSI("T_TRIPLE", name) == approx(value(m.fs.props.temperature_triple[name]))
            # Viscoisty
            assert PropsSI("VISCOSITY", "T", T, "P", P * 1e5, f"PR::{name}") == approx(value(m.fs.state[0].visc_d_phase["Liq"]))
            # Compression Factor
            assert PropsSI("Z", "T", T, "P", P * 1e5, f"PR::{name}") == approx(value(m.fs.state[0].compress_fact_phase["Liq"]))