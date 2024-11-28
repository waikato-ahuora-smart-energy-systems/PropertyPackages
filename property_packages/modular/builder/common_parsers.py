from math import floor
from typing import Any, Dict, List
from .base_parser import BuildBase
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
import csv

class base_units_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return {
            'time': pyunits.s,
            'length': pyunits.m,
            'mass': pyunits.kg,
            'amount': pyunits.mol,
            'temperature': pyunits.K,
        }

class bubble_dew_method_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return LogBubbleDew

class components_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:

        def serialise_component(compound: Compound) -> Dict[str, Any]:
            config = {
                "type": Component,
                "entr_mol_ig_comp": ChemSep,
                "enth_mol_ig_comp": ChemSep,
                "entr_mol_liq_comp": ChemSep,
                "enth_mol_liq_comp": ChemSep,
                "parameter_data": {
                    "mw": (compound["MolecularWeight"].value, pyunits.kg / pyunits.kilomol),
                    "pressure_crit": (compound["CriticalPressure"].value, pyunits.Pa),
                    "temperature_crit": (compound["CriticalTemperature"].value, pyunits.K),
                    "omega": compound["AcentricityFactor"].value,
                }
            }

            keys = [
                ["cp_mol_ig_comp", "RPPHeatCapacityCp"],
                ["pressure_sat_comp", "AntoineVaporPressure"],
                ["dens_mol_liq_comp", "LiquidDensity"],
                ["cp_mol_liq_comp", "LiquidHeatCapacityCp"]
            ]

            for key in keys:
                config[key[0]] = ChemSep
                config["parameter_data"].update({f"{key[0]}_coeff": {
                    "units": compound[key[1]]["units"],
                    "eqno": compound[key[1]]["eqno"],
                    "A": (compound[key[1]]["A"]),
                    "B": (compound[key[1]]["B"]),
                    "C": (compound[key[1]]["C"]),
                    "D": (compound[key[1]]["D"]),
                    "E": (compound[key[1]]["E"]),
                }})

            # Unique logic
            valid_phase = valid_phases(compound)
            if valid_phase != PT.vaporPhase:
                config["phase_equilibrium_form"] = {("Vap", "Liq"): log_fugacity}
            else:
                config["valid_phase_types"] = valid_phase

            config["parameter_data"].update({
                "enth_mol_form_vap_comp_ref": (compound["HeatOfFormation"].value, pyunits.J/pyunits.kilomol),
                "entr_mol_form_vap_comp_ref": (-1 * compound["AbsEntropy"].value, pyunits.J/pyunits.kilomol/pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol / pyunits.K),
            })

            return config

        def valid_phases(compound: Compound) -> PT:
            # Assumption: Anything above hydrogen can exist as both liquid and vapor
            if compound["NormalBoilingPointTemperature"].value >= 21:
                return [PT.liquidPhase, PT.vaporPhase]
            else:
                return PT.vaporPhase
        
        components_output = {}
        for compound in compounds:
            components_output[compound["CompoundID"].value] = serialise_component(compound)
        return components_output

class phase_equilibrium_state_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        if(valid_states.__contains__("Liq") and valid_states.__contains__("Vap")):
            return {("Liq", "Vap"): SmoothVLE}
        else:
            return None

class phases_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
      phases = {}
      for state in valid_states:
          if state == "Liq":
              phases["Liq"] = {
                  "type": LiquidPhase,
                  "equation_of_state": Cubic,
                  "equation_of_state_options": {"type": CubicType.PR},
              }
          elif state == "Vap":
              phases["Vap"] = {
                  "type": VaporPhase,
                  "equation_of_state": Cubic,
                  "equation_of_state_options": {"type": CubicType.PR},
              }
      return phases

class phases_in_equilibrium_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        if(valid_states.__contains__("Liq") and valid_states.__contains__("Vap")):
            return [("Liq", "Vap")]
        else:
            return None

class pressure_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return (101325, pyunits.Pa)

class state_bounds_parser(BuildBase):
    """
    State bounds are used to find optimal solution
    TODO: need to find a way to dynamically determine these
    """
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        min_melting_point = min([compound["NormalMeltingPointTemperature"].value for compound in compounds])
        min_critical_temperature = min([compound["CriticalTemperature"].value for compound in compounds])
        if(min_critical_temperature > 600):
            return {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (min_melting_point, 400, 1000, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            }
        elif(min_critical_temperature > 500):
            return {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (min_melting_point, 300, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            }
        else:
            return {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (min_melting_point, 150, 350, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            }

class state_definition_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> str:
        return FTPx

class temperature_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return (298.15, pyunits.K)

class pr_kappa_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        kappa_parameters = {}
        compound_id_map = {}

        for i, compound1 in enumerate(compounds):
            compound_id_map[str(floor(compound1["LibraryIndex"].value))] = compound1
            for j, compound2 in enumerate(compounds):
                kappa_parameters[(compound1["CompoundID"].value, compound2["CompoundID"].value)] = 0.000
                # Setting all interactions initially to zero

        # TODO: Adjust method so the multiple kappa values for single pair are supported
        # Open and read the interaction data file
        file = open(this_file_dir() + "/data/pr.dat", 'r')
        reader = csv.reader(file, delimiter=';')
        for row in reader:

            # Skip invalid rows
            if len(row) < 4:
                continue
            
            # Extract ID1, ID2, kappa (k12) and comment
            id1, id2, kappa, _ = row

            try:
                kappa_value = float(kappa)
            except ValueError:
                continue

            # Check if the compound IDs exist in the compound list
            if id1 in compound_id_map and id2 in compound_id_map:
                compound1 = compound_id_map[id1]
                compound2 = compound_id_map[id2]
                kappa_parameters[(compound1["CompoundID"].value, compound2["CompoundID"].value)] = kappa_value
                kappa_parameters[(compound2["CompoundID"].value, compound1["CompoundID"].value)] = kappa_value

        return {"PR_kappa": kappa_parameters}
