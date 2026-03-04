import pytest
from ahuora_property_packages.seawater.seawater_builder import build_seawater_package
import pyomo.environ as pe # Pyomo environment
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Heater
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
import matplotlib.pyplot as plt


def test_seawater_heater():
    # Create an empty flowsheet and steam property parameter block.
    model = pe.ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.properties = build_seawater_package(['H2O', 'TDS'])

    # Add a Heater model to the flowsheet.
    model.fs.heater = Heater(property_package=model.fs.properties)

    # Setup the heater model by fixing the inputs and heat duty
    model.fs.heater.inlet.temperature.fix(295)
    model.fs.heater.inlet.pressure.fix(101325)

    sb = model.fs.heater.control_volume.properties_in[0]

    sb.constrain_component(
    sb.mole_frac_comp['H2O']
    , 0.99)

    sb.constrain_component(
    sb.mole_frac_comp['TDS']
    , 0.01)

    sb.constrain_component(
    sb.flow_mol
    , 1)

    model.fs.heater.outlet.temperature.fix(450)

    assert degrees_of_freedom(model) == 0

    model.fs.heater.initialize()

    assert degrees_of_freedom(model) == 0

    # Touch the properties to make sure they are calculated and added to the graph 
    # (since IDAES uses lazy evaluation for properties, they won't be added to the graph until they are accessed)
    sb.flow_mol_phase_comp
    sb.flow_mol
    sb.mole_frac_comp

    solver = pe.SolverFactory("ipopt")
    solver.solve(model,tee=True)

    sb_out = model.fs.heater.control_volume.properties_out[0]

    assert pe.value(sb_out.temperature) == 450
    assert pe.value(model.fs.heater.heat_duty[0]) == pytest.approx(11228.28, rel=1e-4)

