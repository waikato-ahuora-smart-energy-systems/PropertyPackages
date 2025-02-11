# Build and solve a state block.
from property_packages.build_package import build_package
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
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.examples.BT_PR import configuration
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_property import ModularPropertiesInitializer


solver = get_solver(solver="ipopt")

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

def get_m():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)
    iscale.calculate_scaling_factors(m.fs.state[1])
    # applying scaling seems to create convergence issues
    return m

def get_m2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**configuration)
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)
    iscale.calculate_scaling_factors(m.fs.state[1])
    return m

def test_metadata():

    for p in range (50000, 1000000, 50000):
        m = get_m()
        m2 = get_m2()
        sb = m.fs.state[1]
        sb2 = m2.fs.state[1]

        assert_units_consistent(m)
        assert_units_consistent(m2)
        
        sb.flow_mol.fix(1) # What does flow_mol do?
        sb.mole_frac_comp["benzene"].fix(0.5)
        sb.mole_frac_comp["toluene"].fix(0.5)
        sb.pressure.fix(p)
        sb.temperature.fix(300)

        sb2.flow_mol.fix(1)
        sb2.mole_frac_comp["benzene"].fix(0.5)
        sb2.mole_frac_comp["toluene"].fix(0.5)
        sb2.pressure.fix(p)
        sb2.temperature.fix(300)

        m.fs.state.initialize(outlvl=SOUT)
        m2.fs.state.initialize(outlvl=SOUT)

        assert value(sb.flow_vol) == approx(value(sb2.flow_vol), rel=1e-2)
        assert value(sb.flow_mass) == approx(value(sb2.flow_mass), rel=1e-2)
        assert value(sb.phase_frac["Liq"]) == approx(value(sb2.phase_frac["Liq"]), rel=1e-2)
        assert value(sb.phase_frac["Vap"]) == approx(value(sb2.phase_frac["Vap"]), rel=1e-2)

        # assert value(sb._teq["Vap"]) == approx(value(sb2._teq["Vap"]), rel=1e-2)
        assert value(sb.flow_mass_phase["Liq"]) == approx(value(sb2.flow_mass_phase["Liq"]), rel=1e-2)
        assert value(sb.flow_mass_phase["Vap"]) == approx(value(sb2.flow_mass_phase["Vap"]), rel=1e-2)
        assert value(sb.flow_vol_phase["Liq"]) == approx(value(sb2.flow_vol_phase["Liq"]), rel=1e-2)
        assert value(sb.flow_vol_phase["Vap"]) == approx(value(sb2.flow_vol_phase["Vap"]), rel=1e-2)
        assert value(sb.flow_mol_phase["Liq"]) == approx(value(sb2.flow_mol_phase["Liq"]), rel=1e-2)
        assert value(sb.flow_mol_phase["Vap"]) == approx(value(sb2.flow_mol_phase["Vap"]), rel=1e-2)
        assert value(sb.flow_mass_comp["benzene"]) == approx(value(sb2.flow_mass_comp["benzene"]), rel=1e-2)
        assert value(sb.flow_mass_comp["toluene"]) == approx(value(sb2.flow_mass_comp["toluene"]), rel=1e-2)
        assert value(sb.flow_mol_comp["benzene"]) == approx(value(sb2.flow_mol_comp["benzene"]), rel=1e-2)
        assert value(sb.flow_mol_comp["toluene"]) == approx(value(sb2.flow_mol_comp["toluene"]), rel=1e-2)
        assert value(sb.flow_mass_phase_comp["Liq", "benzene"]) == approx(value(sb2.flow_mass_phase_comp["Liq", "benzene"]), rel=1e-2)
        assert value(sb.flow_mass_phase_comp["Vap", "benzene"]) == approx(value(sb2.flow_mass_phase_comp["Vap", "benzene"]), rel=1e-2)
        assert value(sb.flow_mass_phase_comp["Liq", "toluene"]) == approx(value(sb2.flow_mass_phase_comp["Liq", "toluene"]), rel=1e-2)
        assert value(sb.flow_mass_phase_comp["Vap", "toluene"]) == approx(value(sb2.flow_mass_phase_comp["Vap", "toluene"]), rel=1e-2)
        assert value(sb.flow_mol_phase_comp["Liq", "benzene"]) == approx(value(sb2.flow_mol_phase_comp["Liq", "benzene"]), rel=1e-2)
        assert value(sb.flow_mol_phase_comp["Vap", "benzene"]) == approx(value(sb2.flow_mol_phase_comp["Vap", "benzene"]), rel=1e-2)
        assert value(sb.flow_mol_phase_comp["Liq", "toluene"]) == approx(value(sb2.flow_mol_phase_comp["Liq", "toluene"]), rel=1e-2)
        assert value(sb.flow_mol_phase_comp["Vap", "toluene"]) == approx(value(sb2.flow_mol_phase_comp["Vap", "toluene"]), rel=1e-2)
        assert value(sb.mole_frac_comp["benzene"]) == approx(value(sb2.mole_frac_comp["benzene"]), rel=1e-2)
        assert value(sb.mole_frac_comp["toluene"]) == approx(value(sb2.mole_frac_comp["toluene"]), rel=1e-2)

        # Phase compressibility factors (Liq and Vap phases)
        assert value(sb.compress_fact_phase["Liq"]) == approx(value(sb2.compress_fact_phase["Liq"]), rel=1e-1)
        assert value(sb.compress_fact_phase["Vap"]) == approx(value(sb2.compress_fact_phase["Vap"]), rel=1e-1)

        # Concentration and molar concentrations for components in both phases
        assert value(sb.conc_mol_comp["benzene"]) == approx(value(sb2.conc_mol_comp["benzene"]), rel=1e-2)
        assert value(sb.conc_mol_comp["toluene"]) == approx(value(sb2.conc_mol_comp["toluene"]), rel=1e-2)
        assert value(sb.dens_mol_phase["Liq"]) == approx(value(sb2.dens_mol_phase["Liq"]), rel=1e-2)

        # Enthalpy and entropy comparisons
        assert value(sb.enth_mol) == approx(value(sb2.enth_mol), rel=1e-2)
        assert value(sb.enth_mol_phase["Vap"]) == approx(value(sb2.enth_mol_phase["Vap"]), rel=1e-2)
        assert value(sb.enth_mol_phase["Liq"]) == approx(value(sb2.enth_mol_phase["Liq"]), rel=1e-2)
        assert value(sb.enth_mol_phase_comp["Vap", "benzene"]) == approx(value(sb2.enth_mol_phase_comp["Vap", "benzene"]), rel=1e-2)
        assert value(sb.enth_mol_phase_comp["Liq", "benzene"]) == approx(value(sb2.enth_mol_phase_comp["Liq", "benzene"]), rel=1e-2)
        assert value(sb.enth_mol_phase_comp["Vap", "toluene"]) == approx(value(sb2.enth_mol_phase_comp["Vap", "toluene"]), rel=1e-2)
        assert value(sb.enth_mol_phase_comp["Liq", "toluene"]) == approx(value(sb2.enth_mol_phase_comp["Liq", "toluene"]), rel=1e-2)

        assert value(sb.entr_mol) == approx(value(sb2.entr_mol), rel=1e-2)
        assert value(sb.entr_mol_phase["Vap"]) == approx(value(sb2.entr_mol_phase["Vap"]), rel=1e-2)
        assert value(sb.entr_mol_phase["Liq"]) == approx(value(sb2.entr_mol_phase["Liq"]), rel=1e-2)
        assert value(sb.entr_mol_phase_comp["Vap", "benzene"]) == approx(value(sb2.entr_mol_phase_comp["Vap", "benzene"]), rel=1e-2)
        assert value(sb.entr_mol_phase_comp["Liq", "benzene"]) == approx(value(sb2.entr_mol_phase_comp["Liq", "benzene"]), rel=1e-2)
        assert value(sb.entr_mol_phase_comp["Vap", "toluene"]) == approx(value(sb2.entr_mol_phase_comp["Vap", "toluene"]), rel=1e-2)
        assert value(sb.entr_mol_phase_comp["Liq", "toluene"]) == approx(value(sb2.entr_mol_phase_comp["Liq", "toluene"]), rel=1e-2)

        # Fugacity coefficients and Gibbs energy
        assert value(sb.fug_phase_comp["Vap", "benzene"]) == approx(value(sb2.fug_phase_comp["Vap", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.fug_phase_comp["Liq", "benzene"]) == approx(value(sb2.fug_phase_comp["Liq", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.fug_phase_comp["Vap", "toluene"]) == approx(value(sb2.fug_phase_comp["Vap", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.fug_phase_comp["Liq", "toluene"]) == approx(value(sb2.fug_phase_comp["Liq", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.fug_coeff_phase_comp["Vap", "benzene"]) == approx(value(sb2.fug_coeff_phase_comp["Vap", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.fug_coeff_phase_comp["Liq", "benzene"]) == approx(value(sb2.fug_coeff_phase_comp["Liq", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.fug_coeff_phase_comp["Vap", "toluene"]) == approx(value(sb2.fug_coeff_phase_comp["Vap", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.fug_coeff_phase_comp["Liq", "toluene"]) == approx(value(sb2.fug_coeff_phase_comp["Liq", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.gibbs_mol) == approx(value(sb2.gibbs_mol), rel=1e-2)
        assert value(sb.gibbs_mol_phase["Vap"]) == approx(value(sb2.gibbs_mol_phase["Vap"]), rel=1e-2)
        assert value(sb.gibbs_mol_phase["Liq"]) == approx(value(sb2.gibbs_mol_phase["Liq"]), rel=1e-2)
        assert value(sb.gibbs_mol_phase_comp["Vap", "benzene"]) == approx(value(sb2.gibbs_mol_phase_comp["Vap", "benzene"]), rel=1e-2)
        assert value(sb.gibbs_mol_phase_comp["Liq", "benzene"]) == approx(value(sb2.gibbs_mol_phase_comp["Liq", "benzene"]), rel=1e-2)
        assert value(sb.gibbs_mol_phase_comp["Vap", "toluene"]) == approx(value(sb2.gibbs_mol_phase_comp["Vap", "toluene"]), rel=1e-2)
        assert value(sb.gibbs_mol_phase_comp["Liq", "toluene"]) == approx(value(sb2.gibbs_mol_phase_comp["Liq", "toluene"]), rel=1e-2)

        # Mass fractions in phases
        assert value(sb.mass_frac_phase_comp["Vap", "benzene"]) == approx(value(sb2.mass_frac_phase_comp["Vap", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.mass_frac_phase_comp["Liq", "benzene"]) == approx(value(sb2.mass_frac_phase_comp["Liq", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.mass_frac_phase_comp["Vap", "toluene"]) == approx(value(sb2.mass_frac_phase_comp["Vap", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.mass_frac_phase_comp["Liq", "toluene"]) == approx(value(sb2.mass_frac_phase_comp["Liq", "toluene"]), rel=1e-1) # Investigate

        assert value(sb.mw_comp["benzene"]) == approx(value(sb2.mw_comp["benzene"]), rel=1e-2)
        assert value(sb.mw_comp["toluene"]) == approx(value(sb2.mw_comp["toluene"]), rel=1e-2)
        assert value(sb.mw_phase["Vap"]) == approx(value(sb2.mw_phase["Vap"]), rel=1e-2)
        assert value(sb.mw_phase["Liq"]) == approx(value(sb2.mw_phase["Liq"]), rel=1e-2)

        # Other critical pressures and phase pressures
        assert value(sb.pressure_phase_comp["Vap", "benzene"]) == approx(value(sb2.pressure_phase_comp["Vap", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.pressure_phase_comp["Liq", "benzene"]) == approx(value(sb2.pressure_phase_comp["Liq", "benzene"]), rel=1e-1) # Investigate
        assert value(sb.pressure_phase_comp["Vap", "toluene"]) == approx(value(sb2.pressure_phase_comp["Vap", "toluene"]), rel=1e-1) # Investigate
        assert value(sb.pressure_phase_comp["Liq", "toluene"]) == approx(value(sb2.pressure_phase_comp["Liq", "toluene"]), rel=1e-1) # Investigate