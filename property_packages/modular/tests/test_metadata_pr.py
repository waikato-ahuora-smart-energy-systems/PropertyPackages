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

solver = get_solver(solver="ipopt")

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

def get_m():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[1], )
    m.fs.props = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
    iscale.calculate_scaling_factors(m.fs.props)
    iscale.calculate_scaling_factors(m.fs.state[1])
    return m

def test_metadata():

    m = get_m()
    sb = m.fs.state[1]

    assert_units_consistent(m)
    
    sb.flow_mol.fix(100)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    sb.pressure.fix(101325)
    sb.temperature.fix(300)

    # what does init do?
    m.fs.state.initialize(outlvl=SOUT)
    # solver.solve(m, tee=True)

    assert value(sb.flow_vol) == approx(0.009747, rel=1e-3) # Need to verify
    assert value(sb.flow_mass) == approx(8.5125, rel=1e-3) # Need to verify

    # where are these calculated

    assert value(sb.phase_frac["Liq"]) == approx(1.0, rel=1e-3) # Need to verify
    assert value(sb.phase_frac["Vap"]) == approx(0.0, abs=1e-3) # Need to verify

    # assert value(sb._teq["Vap"]) == approx(300, rel=1e-3) # Need to verify

    assert value(sb.flow_mass_phase["Liq"]) == approx(8.512, rel=1e-3)
    assert value(sb.flow_mass_phase["Vap"]) == approx(0.000, abs=1e-3)

    assert value(sb.flow_vol_phase["Liq"]) == approx(0.00974, rel=1e-3)
    assert value(sb.flow_vol_phase["Vap"]) == approx(0.0, abs=1e-3)

    assert value(sb.flow_mol_phase["Liq"]) == approx(100, rel=1e-2)
    assert value(sb.flow_mol_phase["Vap"]) == approx(0.0, abs=1e-2)
    
    assert value(sb.flow_mass_comp["benzene"]) == approx(3.905, rel=1e-3)
    assert value(sb.flow_mass_comp["toluene"]) == approx(4.606, rel=1e-3)

    assert value(sb.flow_mol_comp["benzene"]) == 50
    assert value(sb.flow_mol_comp["toluene"]) == 50

    assert value(sb.flow_mass_phase_comp["Liq", "benzene"]) == approx(3.905, rel=1e-3)
    assert value(sb.flow_mass_phase_comp["Vap", "benzene"]) == approx(0.0, abs=1e-3)
    assert value(sb.flow_mass_phase_comp["Liq", "toluene"]) == approx(4.606, rel=1e-3)
    assert value(sb.flow_mass_phase_comp["Vap", "toluene"]) == approx(0.0, abs=1e-3)

    assert value(sb.flow_mol_phase_comp["Liq", "benzene"]) == approx(50, rel=1e-3)
    assert value(sb.flow_mol_phase_comp["Vap", "benzene"]) == approx(0.0, abs=1e-2)
    assert value(sb.flow_mol_phase_comp["Liq", "toluene"]) == approx(50, rel=1e-3)
    assert value(sb.flow_mol_phase_comp["Vap", "toluene"]) == approx(0.0, abs=1e-2)

    assert value(sb.mole_frac_comp["benzene"]) == 0.5
    assert value(sb.mole_frac_comp["toluene"]) == 0.5

    assert sb.phase_frac is not None
    assert sb.temperature is not None
    assert sb.pressure is not None

    # Not implemented <- check
    # assert value(sb.act_phase_comp) == 1.5
    # assert value(sb.act_phase_comp_true) == 1.3
    # assert value(sb.act_phase_comp_apparent) == 1.4
    # assert value(sb.act_coeff_phase_comp) == 2.5
    # assert value(sb.act_coeff_phase_comp_true) == 2.3
    # assert value(sb.act_coeff_phase_comp_apparent) == 2.4

    assert value(sb.compress_fact_phase["Liq"]) == approx(0.00395, abs=1e-5)
    assert value(sb.compress_fact_phase["Vap"]) == approx(0.949, rel=1e-3)

    # unintialized
    # assert value(sb.compress_fact_crit) == 1.1

    assert value(sb.conc_mol_comp["benzene"]) == approx(5129, rel=1e-3)
    assert value(sb.conc_mol_comp["toluene"]) == approx(5129, rel=1e-3)

    assert value(sb.conc_mol_phase_comp["Liq", "benzene"]) == approx(5129, rel=1e-3)
    assert value(sb.conc_mol_phase_comp["Vap", "benzene"]) == approx(30.3, rel=1e-3) # strange?
    assert value(sb.conc_mol_phase_comp["Liq", "toluene"]) == approx(5129, rel=1e-3)
    assert value(sb.conc_mol_phase_comp["Vap", "toluene"]) == approx(12.479, rel=1e-3) # strange?

    # assert value(sb.conc_mol_phase_comp_apparent) == 0.47 tf is this?
    # assert value(sb.conc_mol_phase_comp_true) == 0.46 recursive aswell?

    # assert value(sb.cp_mass_phase) == 4.5

    # assert value(sb.cp_mol) == 8.5
    # assert value(sb.cp_mol_phase) == 9.0
    # assert value(sb.cp_mol_phase_comp) == 7.5

    # assert value(sb.cv_mass_phase) == 3.5
    # assert value(sb.cv_mol) == 6.5
    # assert value(sb.cv_mol_phase) == 7.0
    # assert value(sb.cv_mol_phase_comp) == 5.5

    # assert value(sb.diffus_phase_comp) == 0.01
    # assert value(sb.diffus_phase_comp["Vap", "benzene"]) == 0.01
    # assert value(sb.diffus_phase_comp["Liq", "toluene"]) == 0.01
    # assert value(sb.diffus_phase_comp["Vap", "toluene"]) == 0.01

    # assert value(sb.diffus_phase_comp_apparent) == 0.015
    # assert value(sb.diffus_phase_comp_true) == 0.02
    # assert value(sb.heat_capacity_ratio_phase) == 1.4
    # assert value(sb.dens_mass) == 800
    # assert value(sb.dens_mass_phase) == 850
    # assert value(sb.dens_mol) == 32
    # assert value(sb.dens_mol_crit) == 29
    # assert value(sb.dens_mol_phase) == 31
    # assert value(sb.energy_internal_mol) == 1000
    # assert value(sb.energy_internal_mol_phase) == 1050
    # assert value(sb.energy_internal_mol_phase_comp) == 1100
    # assert value(sb.enth_mol) == 2000
    # assert value(sb.enth_mol_phase) == 2100
    # assert value(sb.enth_mol_phase_comp) == 2200
    # assert value(sb.entr_mol) == 5.0
    # assert value(sb.entr_mol_phase) == 5.5
    # assert value(sb.entr_mol_phase_comp) == 6.0
    # assert value(sb.fug_phase_comp) == 0.95
    # assert value(sb.fug_coeff_phase_comp) == 0.9
    # assert value(sb.gibbs_mol) == -200
    # assert value(sb.gibbs_mol_phase) == -150
    # assert value(sb.gibbs_mol_phase_comp) == -100
    # assert value(sb.isentropic_speed_sound_phase) == 340
    # assert value(sb.isothermal_speed_sound_phase) == 320
    # assert value(sb.henry) == 0.5
    # assert value(sb.mass_frac_phase_comp) == 0.4
    # assert value(sb.mass_frac_phase_comp_apparent) == 0.42
    # assert value(sb.mass_frac_phase_comp_true) == 0.38
    # assert value(sb.molality_phase_comp) == 1.2
    # assert value(sb.molality_phase_comp_apparent) == 1.3
    # assert value(sb.molality_phase_comp_true) == 1.1
    # assert value(sb.mw) == 100
    # assert value(sb.mw_comp) == 110
    # assert value(sb.mw_phase) == 105
    # assert value(sb.prandtl_number_phase) == 0.7
    # assert value(sb.pressure_crit) == 5000000
    # assert value(sb.pressure_phase_comp) == 101325
    # assert value(sb.pressure_phase_comp_true) == 102000
    # assert value(sb.pressure_phase_comp_apparent) == 101000
    # assert value(sb.pressure_bubble) == 101325
    # assert value(sb.pressure_dew) == 101300
    # assert value(sb.pressure_osm_phase) == 100000
    # assert value(sb.pressure_sat_comp) == 101325
    # assert value(sb.surf_tens_phase) == 0.03
    # assert value(sb.temperature_crit) == 647.1
    # assert value(sb.temperature_bubble) == 373.15
    # assert value(sb.temperature_dew) == 373.65
    # assert value(sb.therm_cond_phase) == 0.02
    # assert value(sb.visc_d_phase) == 0.001
    # assert value(sb.vol_mol_phase) == 0.02
    # assert value(sb.vol_mol_phase_comp) == 0.015
    # assert value(sb.dh_rxn) == -1000
    # assert value(sb.log_act_phase_comp) == 2.5
    # assert value(sb.log_act_phase_solvents) == 2.7
    # assert value(sb.log_act_phase_comp_true) == 2.4
    # assert value(sb.log_act_phase_comp_apparent) == 2.6
    # assert value(sb.log_conc_mol_phase_comp) == 1.2
    # assert value(sb.log_conc_mol_phase_comp_true) == 1.3
    # assert value(sb.log_mass_frac_phase_comp) == 0.8
    # assert value(sb.log_mass_frac_phase_comp_apparent) == 0.85
    # assert value(sb.log_mass_frac_phase_comp_true) == 0.75
    # assert value(sb.log_molality_phase_comp) == 1.1
    # assert value(sb.log_molality_phase_comp_apparent) == 1.2
    # assert value(sb.log_molality_phase_comp_true) == 1.0
    # assert value(sb.log_mole_frac_comp) == 0.5
    # assert value(sb.log_mole_frac_tbub) == 0.6
    # assert value(sb.log_mole_frac_tdew) == 0.55
    # assert value(sb.log_mole_frac_pbub) == 0.7
    # assert value(sb.log_mole_frac_pdew) == 0.65
    # assert value(sb.log_mole_frac_phase_comp) == 0.45
    # assert value(sb.log_mole_frac_phase_comp_apparent) == 0.48
    # assert value(sb.log_mole_frac_phase_comp_true) == 0.42
    # assert value(sb.log_pressure_phase_comp) == 4.5
    # assert value(sb.log_pressure_phase_comp_apparent) == 4.6
    # assert value(sb.log_pressure_phase_comp_true) == 4.7
    # assert value(sb.log_k_eq) == 3.5
    # assert value(sb.flow_mol) == 100
    # # """
