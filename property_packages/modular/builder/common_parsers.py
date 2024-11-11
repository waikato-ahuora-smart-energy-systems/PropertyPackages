from math import floor
from pprint import pprint
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
import csv

"""

    TODO: Need to double check all equation numbers 
    to throw errors if mis-matched or not found

"""


class base_units_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            'time': pyunits.s,
            'length': pyunits.m,
            'mass': pyunits.kg,
            'amount': pyunits.mol,
            'temperature': pyunits.K,
        }

class bubble_dew_method_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return LogBubbleDew

class components_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
    
        def serialise_component(compound: Compound) -> Dict[str, Any]:

            # configuration default to all components 
            config = {
                "type": Component,
                "parameter_data": {
                    "mw": (compound["MolecularWeight"].value, pyunits.kg/pyunits.kilomol),
                    "pressure_crit": (compound["CriticalPressure"].value, pyunits.Pa),
                    "temperature_crit": (compound["CriticalTemperature"].value, pyunits.K),
                    "omega": compound["AcentricityFactor"].value,
                }
            }

            # TODO: update this logic
            if compound["CompoundID"].value == "hydrogen" or compound["CompoundID"].value == "methane":
                config["valid_phase_types"] = PT.vaporPhase
            else:
                config["phase_equilibrium_form"] = {("Vap", "Liq"): log_fugacity}

            # Energies of Formation
            if compound["HeatOfFormation"] is not None: # this does not work when passed
                config["parameter_data"].update({
                    "enth_mol_form_vap_comp_ref": (compound["HeatOfFormation"].value, pyunits.J/pyunits.kilomol)
                })
            else:
                raise ValueError("No Heat of Formation Data")

            if compound["AbsEntropy"] is not None:
                config["parameter_data"].update({
                    "entr_mol_form_vap_comp_ref": (-1 * compound["AbsEntropy"].value, pyunits.J/pyunits.kilomol/pyunits.K)
                })
            else:
                raise ValueError("No Absolute Entropy Data")

            # Ideal Gas Molar Calculations
            # All three properties intrinsically linked together
            if compound["RPPHeatCapacityCp"] is not None:
                if compound["RPPHeatCapacityCp"]["eqno"] == 4:
                    config["enth_mol_ig_comp"] = RPP4
                    config["entr_mol_ig_comp"] = RPP4
                    config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                        "A": (float(compound["RPPHeatCapacityCp"]["A"]), pyunits.J / pyunits.kilomol / pyunits.K),
                        "B": (float(compound["RPPHeatCapacityCp"]["B"]), pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (float(compound["RPPHeatCapacityCp"]["C"]), pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (float(compound["RPPHeatCapacityCp"]["D"]), pyunits.J / pyunits.kilomol / pyunits.K**4),
                    }})
                elif compound["RPPHeatCapacityCp"]["eqno"] == 100 or compound["RPPHeatCapacityCp"]["eqno"] == 5:
                    config["enth_mol_ig_comp"] = ChemSep
                    config["entr_mol_ig_comp"] = ChemSep
                    config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                        "A": (compound["RPPHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K),
                        "B": (compound["RPPHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (compound["RPPHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (compound["RPPHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "E": (compound["RPPHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})
                else:
                    raise ValueError(f"Invalid equation number for heat capacity {compound['RPPHeatCapacityCp']['eqno']}")
            else:
                raise ValueError("No Heat Capacity Data")

            # Saturation Pressure (Vapor)
            if compound["AntoineVaporPressure"] is not None:
                if compound["AntoineVaporPressure"]["eqno"] == 10:
                    config["pressure_sat_comp"] = ChemSep
                    config["parameter_data"].update({"pressure_sat_comp_coeff": {
                        "A": (compound["AntoineVaporPressure"]["A"], None),
                        "B": (compound["AntoineVaporPressure"]["B"], pyunits.K),
                        "C": (compound["AntoineVaporPressure"]["C"], pyunits.K),
                    }})
                else:
                    raise ValueError("No Antoine Vapor Pressure Equation Data")
            else:
                raise ValueError("No Antoine Vapor Pressure Data")

            # Liquid Density
            if compound["LiquidDensity"] is not None:
                if compound["LiquidDensity"]["eqno"] == 105:
                    config["dens_mol_liq_comp"] = Perrys
                    config["parameter_data"].update({"dens_mol_liq_comp_coeff": {
                        "eqn_type": 1,
                        "1": (compound["LiquidDensity"]["A"], pyunits.kmol / pyunits.m**3),
                        "2": (compound["LiquidDensity"]["B"], None),
                        "3": (compound["LiquidDensity"]["C"], pyunits.K),
                        "4": (compound["LiquidDensity"]["D"], None),
                    }})
                else:
                    raise ValueError("No Liquid Density Equation Data")
            else:
                raise ValueError("No Liquid Density Data")

            # Liquid Heat Capacity & Entropy / Enthalpy
            if compound["LiquidHeatCapacityCp"] is not None:
                if compound["LiquidHeatCapacityCp"]["eqno"] == 100:
                    # Uses correct equations to calculate
                    config["enth_mol_liq_comp"] = Perrys
                    config["entr_mol_liq_comp"] = Perrys
                    config["parameter_data"].update({"cp_mol_liq_comp_coeff": {
                        "1": (compound["LiquidHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K**1),
                        "2": (compound["LiquidHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "3": (compound["LiquidHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "4": (compound["LiquidHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "5": (compound["LiquidHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})
                    
                    # ASSUMPTION: Molar heat of formation, liq is zero - given the semi-okay by Ben
                    config["parameter_data"].update({
                        "enth_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol)
                    })

                    config["parameter_data"].update({
                        "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol / pyunits.K)
                    }) 
                elif compound["LiquidHeatCapacityCp"]["eqno"] == 16:
                    # Uses correct equations to calculate
                    config["enth_mol_liq_comp"] = ChemSep
                    config["entr_mol_liq_comp"] = ChemSep
                    config["parameter_data"].update({"cp_mol_liq_comp_coeff": {
                        "A": (compound["LiquidHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K**1),
                        "B": (compound["LiquidHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (compound["LiquidHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (compound["LiquidHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "E": (compound["LiquidHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})

                    # ASSUMPTION: Molar heat of formation, liq is zero - given the semi-okay by Ben
                    config["parameter_data"].update({
                        "enth_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol)
                    })

                    config["parameter_data"].update({
                        "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol / pyunits.K)
                    }) 
                else:
                    raise ValueError(f"No Liquid Heat Capacity Equation Data {compound['LiquidHeatCapacityCp']['eqno']}")
            else: 
                # Compound only exists in vapor phase
                config["parameter_data"].update({"valid_phase_types": PT.vaporphase})

            return config
        
        components_output = {}
        for compound in compounds:
            components_output[compound["CompoundID"].value] = serialise_component(compound)
        return components_output

class phase_equilibrium_state_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {("Vap", "Liq"): SmoothVLE} # TODO: Update to Smooth VLE_v2 https://github.com/IDAES/idaes-pse/blob/main/idaes/models/properties/modular_properties/phase_equil/smooth_VLE_2.py

class phases_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            "Liq": {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
        }

class phases_in_equilibrium_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return [("Vap", "Liq")]

class pressure_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return (101325, pyunits.Pa)

class state_bounds_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (10, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        }

class state_definition_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> str:
        return FTPx

class temperature_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return (298.15, pyunits.K)

class pr_kappa_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        kappa_parameters = {}
        compound_id_map = {}

        for i, compound1 in enumerate(compounds):
            compound_id_map[str(floor(compound1["LibraryIndex"].value))] = compound1
            for j, compound2 in enumerate(compounds):
                kappa_parameters[(compound1["CompoundID"].value, compound2["CompoundID"].value)] = 0.000
                # Setting all interactions initially to zero

        # TODO: Adjust method so the multiple kappa values for single pair are supported
        # TODO: Pass over the correct values
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
                # Convert kappa to a float (if possible)
                kappa_value = float(kappa)
            
            except ValueError:
                continue  # Skip rows with invalid kappa values

            # Check if the compound IDs exist in the compound list

            if id1 in compound_id_map and id2 in compound_id_map:

                compound1 = compound_id_map[id1]
                compound2 = compound_id_map[id2]
                kappa_parameters[(compound1["CompoundID"].value, compound2["CompoundID"].value)] = kappa_value
                kappa_parameters[(compound2["CompoundID"].value, compound1["CompoundID"].value)] = kappa_value

        return {"PR_kappa": kappa_parameters}
