import pyomo.environ as pe# Pyomo environment
from pyomo.environ import value
from idaes.core import FlowsheetBlock, StateBlock
from idaes.models.unit_models import HeatExchanger
from idaes.models.unit_models.heat_exchanger import HX0DInitializer
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback,delta_temperature_amtd_callback
from idaes.models.properties import iapws95
from property_packages.build_package import build_package   
import idaes.logger as idaeslog
from pytest import approx

def test_propane():
    # Tests do not validate data, just used to verify initialisation

    # Create an empty flowsheet and steam property parameter block.
    model = pe.ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.properties = iapws95.Iapws95ParameterBlock()
    model.fs.milk_properties = build_package("milk", ["water", "milk_solid"], ["Liq", "Vap"])

    # Add a Heater model to the flowsheet.
    model.fs.heat_exchanger = HeatExchanger(
        delta_temperature_callback=delta_temperature_lmtd_callback,
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": model.fs.properties},
        tube={"property_package": model.fs.properties}
    )

    model.fs.heat_exchanger.area.fix(1.7)
    model.fs.heat_exchanger.overall_heat_transfer_coefficient[0].fix(1000)
    model.fs.heat_exchanger.shell_inlet.flow_mol.fix(1)
    h = model.fs.properties.htpx(T=(130+273)*pe.units.K, x=1)
    model.fs.heat_exchanger.shell_inlet.pressure.fix(270280)
    model.fs.heat_exchanger.shell_inlet.enth_mol.fix(h)
    model.fs.heat_exchanger.tube_inlet.flow_mol.fix(1)
    model.fs.heat_exchanger.tube_inlet.pressure.fix(101325)

    h1 = model.fs.properties.htpx(T=(12+273)*pe.units.K, p=101325*pe.units.Pa)
    model.fs.heat_exchanger.tube_inlet.enth_mol.fix(h1)
    model.fs.heat_exchanger.initialize(outlvl=idaeslog.INFO_HIGH)

    # Solve the model
    solver = pe.SolverFactory('ipopt')
    results = solver.solve(model, tee=True)

    # Display the results
    model.fs.heat_exchanger.report()
    model.fs.heat_exchanger.cold_side.properties_out[0].display()

    #Unfix area and solve again
    model.fs.heat_exchanger.area.unfix()
    h2 = model.fs.properties.htpx(T=(50+273)*pe.units.K, p=101325*pe.units.Pa)
    model.fs.heat_exchanger.tube_outlet.enth_mol.fix(h2)

    results = solver.solve(model, tee=True)
    model.fs.heat_exchanger.report()
    model.fs.heat_exchanger.cold_side.properties_out[0].display()

def test_milk_custom_props():

    m = pe.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("milk", ["water", "milk_solid"], ["Liq", "Vap"])
    m.fs.state_block = m.fs.properties.build_state_block(defined_state=True)
    n = 0.01
    m.fs.state_block.flow_mol.fix(1)
    m.fs.state_block.temperature.fix(100+273.15)
    m.fs.state_block.pressure.fix(1*100*1000)
    m.fs.state_block.mole_frac_comp["milk_solid"].fix(n)
    m.fs.state_block.mole_frac_comp["water"].fix(1-n)
    m.fs.state_block.initialize()

    m.fs.state_block.enth_mol # trigger enthalpy build
    m.fs.state_block.entr_mol # trigger entropy build

    solver = pe.SolverFactory('ipopt')
    solver.solve(m)

    # Values not verified

    assert value(m.fs.state_block.enth_mol_comp["water"]) == approx(48764, abs=1e3)
    assert value(m.fs.state_block.enth_mol_comp["milk_solid"]) == approx(107, abs=1e3)
    assert value(m.fs.state_block.enth_mol) == approx(48871, abs=1e4)

    assert value(m.fs.state_block.enth_mass_comp["water"]) == approx(2706133, abs=1e5)
    assert value(m.fs.state_block.enth_mass_comp["milk_solid"]) == approx(462, abs=1e4)
    assert value(m.fs.state_block.enth_mass) == approx(2706595, abs=1e5)

    assert value(m.fs.state_block.flow_mol_comp["water"]) == 0.99
    assert value(m.fs.state_block.flow_mol_comp["milk_solid"]) == 0.01

    assert value(m.fs.state_block.flow_mass_comp["water"]) == approx(0.0178, abs=1e-4)
    assert value(m.fs.state_block.flow_mass_comp["milk_solid"]) == approx(0.00232, abs=1e-4)
    assert value(m.fs.state_block.flow_mass) == approx(0.0201, abs=1e-4)

    assert value(m.fs.state_block.entr_mol_comp["water"]) == approx(283, rel=1e-2)
    assert value(m.fs.state_block.entr_mol_comp["milk_solid"]) == approx(196, rel=1e-2)
    assert value(m.fs.state_block.entr_mol) == approx(479, rel=1e-2)

    assert value(m.fs.state_block.entr_mass_comp["water"]) == approx(15712, abs=1e2)
    assert value(m.fs.state_block.entr_mass_comp["milk_solid"]) == approx(847, abs=1)
    assert value(m.fs.state_block.entr_mass) == approx(16560, abs=1e2)

    assert value(m.fs.state_block.vapor_frac) == approx(0.2618, abs=1e-4)
    assert value(m.fs.state_block.flow_vol) == approx(0)
    assert value(m.fs.state_block.total_energy_flow) == approx(54564, abs=1e3)