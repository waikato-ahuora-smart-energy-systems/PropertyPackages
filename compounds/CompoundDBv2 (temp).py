from typing import Dict, Set


class Compound:
    def __init__(self, name: str):
        self.name = name
        self.sources: Dict[str, dict] = {}
        self.supported_packages: Set[str] = set()

    def add_source(self, source_name: str, data: dict):
        self.sources[source_name] = data

    def add_supported_package(self, package_name: str):
        self.supported_packages.add(package_name)

class PropertyPackage:
    def __init__(self, name: str):
        self.name = name
        self.supported_compounds: Set[str] = set()

    def add_supported_compound(self, compound_name: str):
        self.supported_compounds.add(compound_name)

class CompoundRegistry:
    def __init__(self):
        self.compounds: Dict[str, Compound] = {}
        self.packages: Dict[str, PropertyPackage] = {}

    def add_compound_source(self, compound_name: str, source: str, data: dict):
        if compound_name not in self.compounds:
            self.compounds[compound_name] = Compound(compound_name)
        self.compounds[compound_name].add_source(source, data)

    def register_package_support(self, compound_name: str, package_name: str):
        # Ensure both objects exist
        if compound_name not in self.compounds:
            self.compounds[compound_name] = Compound(compound_name)
        if package_name not in self.packages:
            self.packages[package_name] = PropertyPackage(package_name)

        # Register both directions
        self.compounds[compound_name].add_supported_package(package_name)
        self.packages[package_name].add_supported_compound(compound_name)

    def get_supported_packages_for_compound(self, compound_name: str) -> Set[str]:
        return self.compounds[compound_name].supported_packages if compound_name in self.compounds else set()

    def get_supported_compounds_for_package(self, package_name: str) -> Set[str]:
        return self.packages[package_name].supported_compounds if package_name in self.packages else set()


## new workflow for compoundDB

registry = CompoundRegistry()

# Add sources
registry.add_compound_source("Benzene", "genericML", {"MW": 78.11, "BP": 80.1})
registry.add_compound_source("Benzene", "json", {"MW": 78.12, "BP": 80.2})

# Link packages
registry.register_package_support("Benzene", "IDEAL")
registry.register_package_support("Benzene", "PC-SAFT")
registry.register_package_support("Toluene", "IDEAL")

# Query both ways
print(registry.get_supported_packages_for_compound("Benzene"))
print(registry.get_supported_compounds_for_package("IDEAL"))
print(registry.compounds["Benzene"].sources["chemsep"]["MW"])
