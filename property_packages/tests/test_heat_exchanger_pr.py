# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units, Var, Constraint, Expression

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heat_exchanger_ntu import HeatExchangerNTU as HXNTU
from idaes.models.unit_models.heat_exchanger import HeatExchanger, delta_temperature_amtd_callback
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

import idaes.logger as idaeslog

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)


from idaes.models.unit_models.heat_exchanger import (
    delta_temperature_lmtd_callback,
    delta_temperature_lmtd2_callback,
    delta_temperature_amtd_callback,
    delta_temperature_underwood_callback,
    delta_temperature_lmtd_smooth_callback,
    HeatExchanger,
    HeatExchangerFlowPattern,
    HX0DInitializer,
)

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_heat_exchanger_bt():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties1 = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    m.fs.properties2 = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])

    m.fs.unit = HeatExchanger(
        hot_side={"property_package": m.fs.properties1},
        cold_side={"property_package": m.fs.properties2},
        flow_pattern=HeatExchangerFlowPattern.cocurrent,
    )

    m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
    m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
    m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    m.fs.unit.area.fix(1)
    m.fs.unit.overall_heat_transfer_coefficient.fix(100)
    m.fs.unit.cold_side.scaling_factor_pressure = 1

    m.fs.unit.hot_side.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.unit.hot_side.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    m.fs.unit.hot_side.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.unit.hot_side.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    # Important changes, need to set epsilons for new phase equil
    m.fs.unit.cold_side.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.unit.cold_side.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    m.fs.unit.cold_side.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.unit.cold_side.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    m.fs.unit.initialize()

    assert hasattr(m.fs.unit, "hot_side_inlet")
    assert len(m.fs.unit.hot_side_inlet.vars) == 4
    assert hasattr(m.fs.unit.hot_side_inlet, "flow_mol")
    assert hasattr(m.fs.unit.hot_side_inlet, "mole_frac_comp")
    assert hasattr(m.fs.unit.hot_side_inlet, "temperature")
    assert hasattr(m.fs.unit.hot_side_inlet, "pressure")

    assert hasattr(m.fs.unit, "cold_side_inlet")
    assert len(m.fs.unit.cold_side_inlet.vars) == 4
    assert hasattr(m.fs.unit.cold_side_inlet, "flow_mol")
    assert hasattr(m.fs.unit.cold_side_inlet, "mole_frac_comp")
    assert hasattr(m.fs.unit.cold_side_inlet, "temperature")
    assert hasattr(m.fs.unit.cold_side_inlet, "pressure")

    assert hasattr(m.fs.unit, "hot_side_outlet")
    assert len(m.fs.unit.hot_side_outlet.vars) == 4
    assert hasattr(m.fs.unit.hot_side_outlet, "flow_mol")
    assert hasattr(m.fs.unit.hot_side_outlet, "mole_frac_comp")
    assert hasattr(m.fs.unit.hot_side_outlet, "temperature")
    assert hasattr(m.fs.unit.hot_side_outlet, "pressure")

    assert hasattr(m.fs.unit, "cold_side_outlet")
    assert len(m.fs.unit.cold_side_outlet.vars) == 4
    assert hasattr(m.fs.unit.cold_side_outlet, "flow_mol")
    assert hasattr(m.fs.unit.cold_side_outlet, "mole_frac_comp")
    assert hasattr(m.fs.unit.cold_side_outlet, "temperature")
    assert hasattr(m.fs.unit.cold_side_outlet, "pressure")

    assert isinstance(m.fs.unit.overall_heat_transfer_coefficient, Var)
    assert isinstance(m.fs.unit.area, Var)
    assert not hasattr(m.fs.unit, "crossflow_factor")
    assert isinstance(m.fs.unit.heat_duty, Var)
    assert isinstance(m.fs.unit.delta_temperature_in, Var)
    assert isinstance(m.fs.unit.delta_temperature_out, Var)
    assert isinstance(m.fs.unit.unit_heat_balance, Constraint)
    assert isinstance(m.fs.unit.delta_temperature, (Var, Expression))
    assert isinstance(m.fs.unit.heat_transfer_equation, Constraint)

    assert number_variables(m) == 204
    assert number_total_constraints(m) == 90
    assert number_unused_variables(m) == 60

    assert approx(5, abs=1e-3) == value(
        m.fs.unit.hot_side_outlet.flow_mol[0]
    )
    assert approx(359.4, abs=1e-1) == value(
        m.fs.unit.hot_side_outlet.temperature[0]
    )
    assert approx(101325, abs=1e-3) == value(
        m.fs.unit.hot_side_outlet.pressure[0]
    )
    assert approx(1, abs=1e-3) == value(
        m.fs.unit.cold_side_outlet.flow_mol[0]
    )
    assert approx(331.63, abs=1e-1) == value(
        m.fs.unit.cold_side_outlet.temperature[0]
    )
    assert approx(101325, abs=1e-3) == value(
        m.fs.unit.cold_side_outlet.pressure[0]
    )

def test_heat_exchanger_asu():
    pass

    # Need new data

    # m = ConcreteModel()
    # m.fs = FlowsheetBlock(dynamic=False)

    # m.fs.properties1 = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    # m.fs.properties2 = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])

    # m.fs.unit = HeatExchanger(
    #     hot_side={"property_package": m.fs.properties1},
    #     cold_side={"property_package": m.fs.properties2},
    #     flow_pattern=HeatExchangerFlowPattern.cocurrent,
    # )

    # m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
    # m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
    # m.fs.unit.hot_side_inlet.pressure[0].fix(301325)  # Pa
    # m.fs.unit.hot_side_inlet.mole_frac_comp[0, "argon"].fix(0.33)
    # m.fs.unit.hot_side_inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    # m.fs.unit.hot_side_inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)

    # m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
    # m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
    # m.fs.unit.cold_side_inlet.pressure[0].fix(301325)  # Pa
    # m.fs.unit.cold_side_inlet.mole_frac_comp[0, "argon"].fix(0.33)
    # m.fs.unit.cold_side_inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    # m.fs.unit.cold_side_inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)

    # m.fs.unit.area.fix(1)
    # m.fs.unit.overall_heat_transfer_coefficient.fix(100)
    # m.fs.unit.cold_side.scaling_factor_pressure = 1

    # m.fs.unit.hot_side.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    # m.fs.unit.hot_side.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    # m.fs.unit.hot_side.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    # m.fs.unit.hot_side.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    # # Important changes, need to set epsilons for new phase equil
    # m.fs.unit.cold_side.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    # m.fs.unit.cold_side.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    # m.fs.unit.cold_side.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    # m.fs.unit.cold_side.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    # assert degrees_of_freedom(m) == 0

    # m.fs.unit.initialize(outlvl=idaeslog.DEBUG)

    # assert hasattr(m.fs.unit, "hot_side_inlet")
    # assert len(m.fs.unit.hot_side_inlet.vars) == 4
    # assert hasattr(m.fs.unit.hot_side_inlet, "flow_mol")
    # assert hasattr(m.fs.unit.hot_side_inlet, "mole_frac_comp")
    # assert hasattr(m.fs.unit.hot_side_inlet, "temperature")
    # assert hasattr(m.fs.unit.hot_side_inlet, "pressure")

    # assert hasattr(m.fs.unit, "cold_side_inlet")
    # assert len(m.fs.unit.cold_side_inlet.vars) == 4
    # assert hasattr(m.fs.unit.cold_side_inlet, "flow_mol")
    # assert hasattr(m.fs.unit.cold_side_inlet, "mole_frac_comp")
    # assert hasattr(m.fs.unit.cold_side_inlet, "temperature")
    # assert hasattr(m.fs.unit.cold_side_inlet, "pressure")

    # assert hasattr(m.fs.unit, "hot_side_outlet")
    # assert len(m.fs.unit.hot_side_outlet.vars) == 4
    # assert hasattr(m.fs.unit.hot_side_outlet, "flow_mol")
    # assert hasattr(m.fs.unit.hot_side_outlet, "mole_frac_comp")
    # assert hasattr(m.fs.unit.hot_side_outlet, "temperature")
    # assert hasattr(m.fs.unit.hot_side_outlet, "pressure")

    # assert hasattr(m.fs.unit, "cold_side_outlet")
    # assert len(m.fs.unit.cold_side_outlet.vars) == 4
    # assert hasattr(m.fs.unit.cold_side_outlet, "flow_mol")
    # assert hasattr(m.fs.unit.cold_side_outlet, "mole_frac_comp")
    # assert hasattr(m.fs.unit.cold_side_outlet, "temperature")
    # assert hasattr(m.fs.unit.cold_side_outlet, "pressure")

    # assert isinstance(m.fs.unit.overall_heat_transfer_coefficient, Var)
    # assert isinstance(m.fs.unit.area, Var)
    # assert not hasattr(m.fs.unit, "crossflow_factor")
    # assert isinstance(m.fs.unit.heat_duty, Var)
    # assert isinstance(m.fs.unit.delta_temperature_in, Var)
    # assert isinstance(m.fs.unit.delta_temperature_out, Var)
    # assert isinstance(m.fs.unit.unit_heat_balance, Constraint)
    # assert isinstance(m.fs.unit.delta_temperature, (Var, Expression))
    # assert isinstance(m.fs.unit.heat_transfer_equation, Constraint)

    # assert number_variables(m) == 204
    # assert number_total_constraints(m) == 90
    # assert number_unused_variables(m) == 60

    # assert approx(5, abs=1e-3) == value(
    #     m.fs.unit.hot_side_outlet.flow_mol[0]
    # )
    # assert approx(359.4, abs=1e-1) == value(
    #     m.fs.unit.hot_side_outlet.temperature[0]
    # )
    # assert approx(101325, abs=1e-3) == value(
    #     m.fs.unit.hot_side_outlet.pressure[0]
    # )
    # assert approx(1, abs=1e-3) == value(
    #     m.fs.unit.cold_side_outlet.flow_mol[0]
    # )
    # assert approx(331.63, abs=1e-1) == value(
    #     m.fs.unit.cold_side_outlet.temperature[0]
    # )
    # assert approx(101325, abs=1e-3) == value(
    #     m.fs.unit.cold_side_outlet.pressure[0]
    # )