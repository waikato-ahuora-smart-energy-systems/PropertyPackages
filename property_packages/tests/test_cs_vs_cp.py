# Build and solve a state block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
import idaes.core.util.scaling as iscale
from pyomo.environ import assert_optimal_termination

import pytest
from idaes.core.util.model_statistics import (
    degrees_of_freedom)

solver = get_solver(solver="ipopt")

from numpy import logspace
import idaes.logger as idaeslog
SOUT = idaeslog.INFO

from idaes.models.properties.modular_properties.coolprop.coolprop_wrapper import (
    CoolPropWrapper,
    CoolPropExpressionError,
    CoolPropPropertyError,
)
from pyomo.common.dependencies import attempt_import

from compounds.CompoundDB import get_compound_names, get_compound

from tabulate import tabulate

CoolProp, coolprop_available = attempt_import("CoolProp.CoolProp")

# keys = [
#     "flow_mol", "flow_vol", "flow_mass", "flow_mass_phase", "flow_vol_phase", "flow_mol_phase", 
#     "flow_mass_comp", "flow_mol_comp", "flow_mass_phase_comp", "flow_mol_phase_comp", "mole_frac_comp", 
#     "mole_frac_phase_comp", "phase_frac", "temperature", "pressure", "act_phase_comp", 
#     "act_phase_comp_true", "act_phase_comp_apparent", "act_coeff_phase_comp", "act_coeff_phase_comp_true", 
#     "act_coeff_phase_comp_apparent", "compress_fact_phase", "compress_fact_crit", "conc_mol_comp", 
#     "conc_mol_phase_comp", "conc_mol_phase_comp_apparent", "conc_mol_phase_comp_true", "cp_mass_phase", 
#     "cp_mol", "cp_mol_phase", "cp_mol_phase_comp", "cv_mass_phase", "cv_mol", "cv_mol_phase", 
#     "cv_mol_phase_comp", "diffus_phase_comp", "diffus_phase_comp_apparent", "diffus_phase_comp_true", 
#     "heat_capacity_ratio_phase", "dens_mass", "dens_mass_phase", "dens_mol", "dens_mol_crit", "dens_mol_phase", 
#     "energy_internal_mol", "energy_internal_mol_phase", "energy_internal_mol_phase_comp", "enth_mol", 
#     "enth_mol_phase", "enth_mol_phase_comp", "entr_mol", "entr_mol_phase", "entr_mol_phase_comp", 
#     "fug_phase_comp", "fug_coeff_phase_comp", "gibbs_mol", "gibbs_mol_phase", "gibbs_mol_phase_comp", 
#     "isentropic_speed_sound_phase", "isothermal_speed_sound_phase", "henry", "mass_frac_phase_comp", 
#     "mass_frac_phase_comp_apparent", "mass_frac_phase_comp_true", "molality_phase_comp", 
#     "molality_phase_comp_apparent", "molality_phase_comp_true", "mw", "mw_comp", "mw_phase", "prandtl_number"
# ]

def relative_error(a, e):
    if a == 0 or e == 0:
        return "NA"
    return f"{abs((a - e)/(e))*100:.3f}%"
    #return {a, e}



def test_cs_vs_cp():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = build_package("peng-robinson", ["benzene"], ["Liq", "Vap"])
    m.fs.state = m.fs.properties.build_state_block([1], defined_state=True)
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

    iscale.calculate_scaling_factors(m.fs.properties)
    iscale.calculate_scaling_factors(m.fs.state[1])

    m.fs.state[1].flow_mol.fix(1)
    m.fs.state[1].pressure.fix(201325)
    m.fs.state[1].temperature.fix(298)
    m.fs.state[1].mole_frac_comp["benzene"].fix(1)

    assert degrees_of_freedom(m) == 0

    m.fs.state.initialize(outlvl=idaeslog.WARNING)

    assert value(m.fs.state[1].enth_mol) == approx(50290, rel=1e-2)

    for c in get_compound_names():
        try:
            CoolPropWrapper._load_component_data(c)
        except:
            #print(f"Component {c} not found in CoolProp")  
            continue
        else:
            cp_tc = CoolPropWrapper.get_parameter_value(c, "temperature_crit")[0]
            cs_tc = get_compound(c)["CriticalTemperature"].value

            cp_pc = CoolPropWrapper.get_parameter_value(c, "pressure_crit")[0] / 1e5
            cs_pc = get_compound(c)["CriticalPressure"].value / 1e5

            cp_mw = CoolPropWrapper.get_parameter_value(c, "mw")[0] * 1e3
            cs_mw = get_compound(c)["MolecularWeight"].value

            cp_omega = CoolPropWrapper.get_parameter_value(c, "omega")[0]
            cs_omega = get_compound(c)["AcentricityFactor"].value

            print(f"Compound: {c}")

            print(tabulate([
                ["Tc", cs_tc, cp_tc, relative_error(cs_tc, cp_tc)], 
                ["Pc", cs_pc, cp_pc, relative_error(cs_pc, cp_pc)], 
                ["mw", cs_mw, cp_mw, relative_error(cs_mw, cp_mw)],
                ["omega", cs_omega, cp_omega, relative_error(cs_omega, cp_omega)]
            ], ["prop", "chem sep", "coolprop", "rel err %"], tablefmt="grid"))

            # assert CoolPropWrapper.get_parameter_value(c, "temperature_crit")[0] == approx(get_compound(c)["CriticalTemperature"].value, rel=1e-3)
            # assert CoolPropWrapper.get_parameter_value(c, "pressure_crit")[0] == approx(get_compound(c)["CriticalPressure"].value, rel=1e-3)
            # assert CoolPropWrapper.get_parameter_value(c, "mw")[0] * 1e3 == approx(get_compound(c)["MolecularWeight"].value, rel=1e-3)
            # assert CoolPropWrapper.get_parameter_value(c, "omega")[0] == approx(get_compound(c)["AcentricityFactor"].value, rel=1e-2)