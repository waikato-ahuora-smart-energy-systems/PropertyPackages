from pprint import pprint
import io



class Generator:
    """

    NOTE: units and other marked up text elements should avoid using spaces
    as pprint may interpret them as newlines.

    """

    def __init__(self, reader):
        self.compounds = reader.get_all_compounds()
        self.reader = reader
        self.phases = ["Vap", "Liq"] # Hardcoded for now
        self.configuration = {}

    def _build_configuration(self):

        self.configuration.update({"components": self._build_component_configuration()})
        self.configuration.update({"parameter_data": self._build_interaction_parameters()})

        # Default Configuration
        self.configuration.update({"phases": self._build_phases()})
        self.configuration.update({"base_units": self._build_base_units()})
        self.configuration.update({"state_bounds": self._build_state_bounds()})
        self.configuration.update({"state_definition": "<text>FTPx</text>"})
        self.configuration.update({"pressure_ref": (101325, "<text>pyunits.Pa</text>")})
        self.configuration.update({"temperature_ref": (298.15, "<text>pyunits.K</text>")})
        self.configuration.update({"phases_in_equilibrium": [("Vap", "Liq")]})
        self.configuration.update({"phase_equilibrium_state": {("Vap", "Liq"): "<text>SmoothVLE</text>"}})
        self.configuration.update({"bubble_dew_method": "<text>LogBubbleDew</text>"})

    def _output_configuration(self):
        # Path to the output file
        file_path = 'config.py'

        output = io.StringIO()
        pprint(self.configuration, stream=output)
        output_str = output.getvalue()
        prefixed_output = f'''
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.pure import RPP4, Perrys
from idaes.core import LiquidPhase, VaporPhase, Component
from pyomo.environ import units as pyunits
from mod.chem_sep import ChemSep

config = {output_str}
        '''

        # Removing markup language
        cleaned = prefixed_output.replace("'<text>", "")
        cleaned = cleaned.replace("</text>'", "")
        with open(file_path, 'w') as file2:
            file2.write(cleaned)

        print("Markup language removed")

        print(f"Dictionary has been written to {file_path}")


    def get_configuration(self):
        if self.configuration:
            return self.configuration
        else:
            self._build_configuration()
            return self.configuration

    def _build_component_configuration(self):
        """
        - Returns a transmissible property package configuration dictionary
        """
        components = {}

        for compound in self.compounds:
            if self.reader.does_param_exist(compound):
                # Assuming these must exist, are formatted the
                # same or are constant for all compounds
                config = {
                    "type": "<text>Component</text>",
                    "phase_equilibrium_form": {("Vap", "Liq"): "<text>log_fugacity</text>"},
                    "parameter_data": {
                        "mw": (self.reader.get_property(compound, "MolecularWeight"), "<text>pyunits.kg/pyunits.kilomol</text>"),
                        "pressure_crit": (self.reader.get_property(compound, "CriticalPressure"), "<text>pyunits.Pa</text>"),
                        "temperature_crit": (self.reader.get_property(compound, "CriticalTemperature"), "<text>pyunits.K</text>"),
                        "omega": self.reader.get_property(compound, "AcentricityFactor"),
                    }
                }

                # Energies of Formation
                if self.reader.get_property(compound, "HeatOfFormation") is not None:
                    config["parameter_data"].update({
                        "enth_mol_form_vap_comp_ref": (float(self.reader.get_property(compound, "HeatOfFormation")), "<text>pyunits.J/pyunits.kilomol</text>")
                    })
                if self.reader.get_property(compound, "AbsEntropy") is not None:
                    config["parameter_data"].update({
                        "entr_mol_form_vap_comp_ref": (-1 * self.reader.get_property(compound, "AbsEntropy"), "<text>pyunits.J/pyunits.kilomol/pyunits.K</text>")
                    })

                # Ideal Gas Molar Calculations
                # All three properties intrinsically linked together
                if not (self.reader.get_coeff(compound, "RPPHeatCapacityCp") is None):
                    if int(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["eqno"]) == 4:
                        config["enth_mol_ig_comp"] = "<text>RPP4</text>"
                        config["entr_mol_ig_comp"] = "<text>RPP4</text>"
                        config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                            "A": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["A"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K</text>"),
                            "B": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["B"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**2</text>"),
                            "C": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["C"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**3</text>"),
                            "D": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["D"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**4</text>"),
                        }})
                    elif int(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["eqno"]) == 100:
                        config["enth_mol_ig_comp"] = "<text>ChemSep</text>"
                        config["entr_mol_ig_comp"] = "<text>ChemSep</text>"
                        config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                            "A": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["A"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K</text>"),
                            "B": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["B"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**2</text>"),
                            "C": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["C"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**3</text>"),
                            "D": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["D"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**4</text>"),
                            "E": (float(self.reader.get_coeff(compound, "RPPHeatCapacityCp")["E"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**5</text>"),
                        }})

                # Saturation Pressure (Vapor)
                if self.reader.get_coeff(compound, "VaporPressure") is not None:
                    config["pressure_sat_comp"] = "<text>ChemSep</text>"
                    config["parameter_data"].update({"pressure_sat_comp_coeff": {
                        "A": (float(self.reader.get_coeff(compound, "VaporPressure")["A"]), "<text>None</text>"),
                        "B": (float(self.reader.get_coeff(compound, "VaporPressure")["B"]), "<text>pyunits.K</text>"),
                        "C": (float(self.reader.get_coeff(compound, "VaporPressure")["C"]), "<text>pyunits.K</text>"),
                    }})

                # Liquid Density
                if self.reader.get_coeff(compound, "LiquidDensity") is not None:
                    if int(self.reader.get_coeff(compound, "LiquidDensity")["eqno"]) == 105:
                        config["dens_mol_liq_comp"] = "<text>Perrys</text>"
                        config["parameter_data"].update({"dens_mol_liq_comp_coeff": {
                            "eqn_type": 1,
                            "1": (float(self.reader.get_coeff(compound, "LiquidDensity")["A"]), "<text>pyunits.kmol/pyunits.m**3</text>"),
                            "2": (float(self.reader.get_coeff(compound, "LiquidDensity")["B"]), "<text>None</text>"),
                            "3": (float(self.reader.get_coeff(compound, "LiquidDensity")["C"]), "<text>pyunits.K</text>"),
                            "4": (float(self.reader.get_coeff(compound, "LiquidDensity")["D"]), "<text>None</text>"),
                        }})

                # Liquid Heat Capacity & Entropy / Enthalpy
                if self.reader.get_coeff(compound, "LiquidHeatCapacityCp") is not None:
                    if int(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["eqno"]) == 100:
                        # Uses correct equations to calculate
                        config["enth_mol_liq_comp"] = "<text>Perrys</text>"
                        config["entr_mol_liq_comp"] = "<text>Perrys</text>"
                        config["parameter_data"].update({"cp_mol_liq_comp_coeff": {
                            "1": (float(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["A"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**1</text>"),
                            "2": (float(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["B"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**2</text>"),
                            "3": (float(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["C"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**3</text>"),
                            "4": (float(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["D"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**4</text>"),
                            "5": (float(self.reader.get_coeff(compound, "LiquidHeatCapacityCp")["E"]), "<text>pyunits.J/pyunits.kilomol/pyunits.K**5</text>"),
                        }})

                        # ASSUMPTION: Molar heat of formation, liq is zero - given the semi-okay by Ben
                        config["parameter_data"].update({
                            "enth_mol_form_liq_comp_ref": (0, "<text>pyunits.J/pyunits.kilomol</text>")
                        })
                        config["parameter_data"].update({
                            "entr_mol_form_liq_comp_ref": (0, "<text>pyunits.J/pyunits.kilomol/pyunits.K</text>")
                        })

                components.update({compound.lower(): config})

        return components

    def _build_phases(self):
        config = {}
        for phase in self.phases:
            if phase.lower() == "vap":
                config.update({"Vap": self._build_vapor_phase()})
            elif phase.lower() == "liq":
                config.update({"Liq": self._build_liquid_phase()})
        return config

    def _build_liquid_phase(self):
        return {
            "type": "<text>LiquidPhase</text>",
            "equation_of_state": "<text>Cubic</text>",
            "equation_of_state_options": {"type": "<text>CubicType.PR</text>"},
        }

    def _build_vapor_phase(self):
        return {
            "type": "<text>VaporPhase</text>",
            "equation_of_state": "<text>Cubic</text>",
            "equation_of_state_options": {"type": "<text>CubicType.PR</text>"},
        }


    def _does_interaction_parameter_exist(self, compound1, compound2):
        """
        Checks if a given interaction parameter exists between
        two compounds.
        """
        if self.reader.get_param(compound1, compound2) is not None:
            return True
        elif self.reader.get_param(compound2, compound1) is not None:
            return True
        else:
            return False

    def _build_interaction_parameters(self):
        """
        Method assumes that if only a single kappa value exists,
        compounds interact "symmetrically" i.e. ratio from C1 : C2
        is equal to C2 : C1. Finally, if compound pair does not exist
        we assume that the kappa value is 0.000.
        """

        config = {
            "PR_kappa": {}
        }

        unique_keys = {}

        for key in self.reader.param_data.keys():

            # key example = '906;902'

            keys = key.split(";")
            value = self.reader.param_data[key]

            compound1 = self.reader.get_compound_name(keys[0])
            compound2 = self.reader.get_compound_name(keys[1])

            if compound1 is None or compound2 is None:
                continue

            compound1 = compound1.lower()
            compound2 = compound2.lower()

            unique_keys[compound1] = True
            unique_keys[compound2] = True

            # getting compound names

            config["PR_kappa"].update({
                (compound1, compound1): 0,
                (compound2, compound2): 0,
                (compound2, compound1): float(value),
                (compound1, compound2): float(value)})

        # setting all other parameters to zero
        for key1 in unique_keys:
            for key2 in unique_keys:
                if (key1, key2) not in config["PR_kappa"]:
                    config["PR_kappa"].update({
                        (key1, key2): 0,
                        (key2, key1): 0,
                    })

        return config

    def _build_base_units(self):
        # Default units are SI Molar
        return {
            "time": "<text>pyunits.s</text>",
            "length": "<text>pyunits.m</text>",
            "mass": "<text>pyunits.kg</text>",
            "temperature": "<text>pyunits.K</text>",
            "amount": "<text>pyunits.mol</text>",
        }

    def _build_state_bounds(self):
        # Default state bounds
        return {
            "flow_mol": (0, 100, 1000, "<text>pyunits.mol / pyunits.s</text>"),
            "temperature": (273.15, 300, 500, "<text>pyunits.K</text>"),
            "pressure": (5e4, 1e5, 1e6, "<text>pyunits.Pa</text>"),
        }
