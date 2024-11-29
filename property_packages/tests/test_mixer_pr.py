# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heater import Heater
from idaes.core.util.tables import _get_state_from_port

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the mixer unit model
from idaes.models.unit_models import Mixer, MomentumMixingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the BTX_ideal property package to create a properties block for the flowsheet
from idaes.models.properties.activity_coeff_models import BTX_activity_coeff_VLE

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import _get_state_from_port

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_mixer():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    
    m.fs.mixer_1 = Mixer(property_package=m.fs.properties,
                        num_inlets=2,
                        momentum_mixing_type=MomentumMixingType.minimize)
    
    assert degrees_of_freedom(m) == 10

    # Benzene stream
    m.fs.mixer_1.inlet_1.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_1.inlet_1.mole_frac_comp[0, "benzene"].fix(0.999)
    m.fs.mixer_1.inlet_1.mole_frac_comp[0, "toluene"].fix(0.001)
    m.fs.mixer_1.inlet_1.pressure.fix(101325*2) # Pa
    m.fs.mixer_1.inlet_1.temperature.fix(353) # K

    # Toluene stream
    m.fs.mixer_1.inlet_2.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_1.inlet_2.mole_frac_comp[0, "benzene"].fix(0.001)
    m.fs.mixer_1.inlet_2.mole_frac_comp[0, "toluene"].fix(0.999)
    m.fs.mixer_1.inlet_2.pressure.fix(101325*4) # Pa
    m.fs.mixer_1.inlet_2.temperature.fix(356) # K

    assert degrees_of_freedom(m) == 0

    m.fs.mixer_1.initialize(outlvl=idaeslog.WARNING)

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    inlet_1 = m.fs.mixer_1.inlet_1_state[0]
    inlet_2 = m.fs.mixer_1.inlet_2_state[0]
    outlet = m.fs.mixer_1.mixed_state[0]

    assert approx(outlet.mole_frac_comp["benzene"].value) == 0.5
    assert approx(outlet.mole_frac_comp["toluene"].value) == 0.5

    assert value(outlet.flow_mol_phase["Vap"]) == approx(0, abs=1e-2)
    assert value(outlet.flow_mol_phase["Liq"]) == approx(200, abs=1e-2)
    assert value(inlet_1.flow_mol_phase["Vap"]) == approx(0, abs=1e-2)
    assert value(inlet_1.flow_mol_phase["Liq"]) == approx(100, abs=1e-2)
    assert value(inlet_2.flow_mol_phase["Vap"]) == approx(0, abs=1e-2)
    assert value(inlet_2.flow_mol_phase["Liq"]) == approx(100, abs=1e-2)
    assert value(outlet.temperature) == approx(354.5, abs=1e-1)
    assert value(outlet.pressure) == approx(101325*2, abs=1e-2)


    m.fs.mixer_1.report()