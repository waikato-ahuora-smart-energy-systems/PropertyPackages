from idaes.core.base.process_base import ProcessBaseBlock
import xml.etree.ElementTree as ET


class XMLReader:

    def load(self):
        """
        Loads raw chemsep xml information into class
        """

        print("Loading ChemSep XML Data...")

        # adding both compound sets together and saving into raw_data
        self.chem_sep_1 = ET.parse(open("./data/chemsep/chemsep1.xml"))
        self.chem_sep_2 = ET.parse(open("./data/chemsep/chemsep2.xml"))
        self.chem_sep_1.getroot().extend(list(self.chem_sep_2.getroot()))
        self.raw_data = self.chem_sep_1.getroot()
        self.available_compounds = {}
        self.param_data = {}

        # adding pr information into pr_data
        self.pr_data = open("./data/chemsep/interaction_params.dat").readlines()[1:] # removing first and last line

        # creating a dictionary for pair values
        for line in self.pr_data:
            line = line.split(";")
            self.param_data[line[0] + ";" + line[1]] = line[2]
            self.available_compounds[line[0]] = True
            self.available_compounds[line[1]] = True

        print("ChemSep XML Data Loaded.")

    def create_configuration(self):
        """
        Creates an executable python script for the idaes
        modular property package framework with the current
        chemsep data and equations
        """

        # looping through all compounds


    def create_compound(self, compound_id):
        """
        Creates a compound object with all the properties
        for a given compound ID
        """

    def get_property(self, compound_name, prop_name):
        """
        Retrieves the value of a specific property for a given compound ID.

        Parameters:
        - compound_id (str): The ID of the compound.
        - prop_name (str): The name of the property.

        Returns:
        - str or None: The value of the property if it exists, None otherwise.
        """

        xml_compound = self.raw_data.find(f'.//CompoundID[@name="Name"][@value="{compound_name}"]...')
        xml_property = xml_compound.find(f'.//{prop_name}')

        # Checking if property exists
        if xml_property is None:
            return None

        try:
            float(xml_property.get("value"))
            return float(xml_property.get("value"))
        except ValueError:
            return xml_property.get("value")

    def get_compound_name(self, compound_id):
        xml_compound1 = self.raw_data.find(f'.//LibraryIndex[@name="Index"][@value="{compound_id}"]...')
        if xml_compound1 is None:
            return None
        return xml_compound1.find(f'.//CompoundID').get("value")

    def get_param(self, compound1, compound2):
        """

        """

        xml_compound1 = self.raw_data.find(f'.//CompoundID[@name="Name"][@value="{compound1}"]...')
        xml_property1 = xml_compound1.find(f'.//LibraryIndex').get("value")
        xml_compound2 = self.raw_data.find(f'.//CompoundID[@name="Name"][@value="{compound2}"]...')
        xml_property2 = xml_compound2.find(f'.//LibraryIndex').get("value")

        # Checking if property exists
        if xml_property1 is None or xml_property2 is None:
            return None

        # we assume that properties are symmetric
        # we assume that if the property doesn't exist
        # assuming that the ratio of the interaction parameter to itself is 0

        if self.param_data[xml_property1 + ";" + xml_property2] is None:
            return None
        else:
            return self.param_data[xml_property1 + ";" + xml_property2]

    def does_param_exist(self, compound):
        xml_compound1 = self.raw_data.find(f'.//CompoundID[@name="Name"][@value="{compound}"]...')
        xml_property1 = xml_compound1.find(f'.//LibraryIndex').get("value")
        return xml_property1 in self.available_compounds

    def get_coeff(self, compound_name, prop_name):
        """
        Retrieves the coefficients for a given compound ID and property name.

        Parameters:
        - compound_id (str): The ID of the compound.
        - prop_name (str): The name of the property.

        Returns:
        - results (dict): A dictionary containing the coefficients. The keys are the properties
          ('A', 'B', 'C', 'D', 'E', 'T_Min', 'T_Max', 'Equation') and the values are the corresponding
          coefficient values.
        """

        xml_compound = self.raw_data.find(f'.//CompoundID[@name="Name"][@value="{compound_name}"]...')
        xml_coefficient = xml_compound.find(f'.//{prop_name}')
        results = {}
        if xml_coefficient is not None:
            arr = ['A', 'B', 'C', 'D', 'E', 'T_Min', 'T_Max', 'eqno']
            for prop in arr:
                obj = xml_coefficient.find(f'.//{prop}')
                if obj is not None:
                    results[prop] = obj.get("value")
        if results == {}:
            return None
        return results

    def get_all_compounds(self):
        """
        Retrieves all compound IDs from the XML data.

        Returns:
        - list: A list of all compound IDs.
        """
        names = self.raw_data.findall('.//CompoundID[@name="Name"]')
        return [compound.get("value") for compound in names]

    def test_compound_name(self):
        self.compare_properties("Name", "CompoundID")

    def test_compound_formula(self):
        self.compare_properties("Formula", "StructureFormula")

    def test_compound_family(self):
        self.compare_properties("ChemSepFamily", "Family")

    def test_compound_critical_temperature(self):
        self.compare_properties("Critical_Temperature", "CriticalTemperature")

    def test_compound_critical_pressure(self):
        self.compare_properties("Critical_Pressure", "CriticalPressure")

    def test_compound_critical_volume(self):
        self.compare_properties("Critical_Volume", "CriticalVolume")

    def test_compound_critical_compressibility(self):
        self.compare_properties("Critical_Compressibility", "CriticalCompressibility")

    def test_compound_NBPT(self):
        self.compare_properties("Normal_Boiling_Point", "NormalBoilingPointTemperature")

    def test_compound_NMPT(self):
        self.compare_properties("NormalMeltingPointTemperature", "NormalMeltingPointTemperature")

    def test_compound_TPT(self):
        self.compare_properties("TriplePointTemperature", "TriplePointTemperature")

    def test_compound_TPP(self):
        self.compare_properties("TriplePointPressure", "TriplePointPressure")

    def test_compound_molecular_weight(self):
        self.compare_properties("Molar_Weight", "MolecularWeight")

    def test_compound_liquid_volume_at_NBP(self):
        self.compare_properties("LiquidVolumeAtNormalBoilingPoint", "LiquidVolumeAtNormalBoilingPoint")

    def test_compound_acentric_factor(self):
        self.compare_properties("Acentric_Factor", "AcentricityFactor")

    def test_compound_radius_of_gyration(self):
        self.compare_properties("RadiusOfGyration", "RadiusOfGyration")

    def test_compound_solubility_parameter(self):
        self.compare_properties("SolubilityParameter", "SolubilityParameter")

    def test_compound_dipole_moment(self):
        self.compare_properties("Dipole_Moment", "DipoleMoment")

    def test_compound_van_der_waals_volume(self):
        self.compare_properties("VanDerWaalsVolume", "VanDerWaalsVolume")

    def test_compound_van_der_waals_area(self):
        self.compare_properties("VanDerWaalsArea", "VanDerWaalsArea")

    def test_compound_heat_of_formation(self):
        self.compare_properties("IG_Enthalpy_of_Formation_25C", "HeatOfFormation")

    def test_compound_IG_energy_of_formation(self):
        self.compare_properties("IG_Entropy_of_Formation_25C", "AbsEntropy")

    def test_compound_IG_absolute_entropy(self):
        self.compare_properties("IG_Gibbs_Energy_of_Formation_25C", "GibbsEnergyOfFormation")

    def test_compound_heat_of_fusion_at_melting_point(self):
        self.compare_properties("TemperatureOfFusion", "HeatOfFusionAtMeltingPoint")

    def test_compound_matthias_copeman_c1(self):
        self.compare_properties("MatthiasCopemanC1", "MatthiasCopemanC1")

    def test_compound_heat_of_combustion(self):
        self.compare_properties("HeatOfCombustion", "HeatOfCombustion")

    def test_compound_solid_density(self):
        self.compare_coefficients("SolidDensity")

    def test_compound_liquid_density(self):
        self.compare_coefficients("LiquidDensity")

    def test_compound_vapor_pressure(self):
        self.compare_coefficients("VaporPressure")

    def test_compound_heat_of_vaporization(self):
        self.compare_coefficients("HeatOfVaporization")

    def test_compound_liquid_heat_capacity_cp(self):
        self.compare_coefficients("LiquidHeatCapacityCp")

    def test_compound_ideal_gas_heat_capacity_cp(self):
        self.compare_coefficients("IdealGasHeatCapacityCp")

    def test_compound_second_virial_coefficient(self):
        self.compare_coefficients("SecondViricalCoefficient")

    def test_compound_liquid_viscosity(self):
        self.compare_coefficients("LiquidViscosity")

    def test_compound_vapor_viscosity(self):
        self.compare_coefficients("VaporViscosity")

    def test_compound_liquid_thermal_conductivity(self):
        self.compare_coefficients("LiquidThermalConductivity")

    def test_compound_vapor_thermal_conductivity(self):
        self.compare_coefficients("VaporThermalConductivity")

    def test_compound_surface_tension(self):
        self.compare_coefficients("SurfaceTension")

    def test_compound_RPP_heat_capacity_cp(self):
        self.compare_coefficients("RPPHeatCapacityCp")

    def test_compound_relative_static_permittivity(self):
        self.compare_coefficients("RelativeStaticPermittivity")

    def test_compound_antoine_vapor_pressure(self):
        self.compare_coefficients("AntoineVaporPressure")

    def test_compound_liquid_viscosity_RPS(self):
        self.compare_coefficients("LiquidViscosityRPS")

    def test_liquid_heat_capacity_regressed(self):
        self.compare_coefficients("LiquidHeatCapacityRegressed")
