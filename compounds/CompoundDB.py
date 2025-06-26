from typing import Dict, Set
from abc import ABC, abstractmethod

__author__ = "Mahaki Leach"


class Compound:
    def __init__(self, name: str):
        self.__name = name
        self.__sources: Dict[str, dict] = {}
        self.__packages: Set[str] = set() # supported packages

    # Read only properties
    @property
    def name(self):
        return self.__name

    @property
    def sources(self):
        return self.__sources

    @property
    def packages(self):
        return self.__packages
    
    def get_source(self, source_name: str) -> dict:
        """
        Get the data from a specific source.
        
        Args:
            source_name (str): Name of the source to retrieve data from.
        
        Returns:
            dict: Data associated with the source if it exists, otherwise None.
        """
        return self.sources.get(source_name, None)
    
    def add_source(self, source_name: str, data: dict):
        if source_name in self.sources:
            raise ValueError(f"Source {source_name} already exists for compound {self.name}.")
        self.sources[source_name] = data

    def add_package(self, package_name: str):
        self.packages.add(package_name)


class PropertyPackage(ABC):

    def __init__(self, name: str):
        # Unique name for the package
        self.__name = name
        # Set of compounds supported by this package
        self.__compounds = set()
    
    # Read only properties
    @property
    def compounds(self):
        return self.__compounds
    
    @property
    def name(self):
        return self.__name

    # Generic methods
    def register_compound(self, name: str):
        self.compounds.add(name) 
    
    def unregister_compound(self, name: str):
        if name in self.compounds:
            self.compounds.remove(name)
        else:
            raise ValueError(f"Compound {name} is not registered in package {self.name}")
    
    def check_registered(self, name: str) -> bool:
        """
        Check if a compound is registered in this package.
        
        Args:
            name (str): Name of the compound to check.
        
        Returns:
            bool: True if the compound is registered, False otherwise.
        """
        return name in self.compounds
    
    def check_supported_compounds(self, compounds: Set[str], strict: bool=True):
        """
        Check if a set of compounds are supported by this package.
        
        Args:
            names (Set[str]): Name of the compound to check.
            strict (bool):  If True, all compounds must be supported.
                            If False, only a single compound must be supported.
        
        Returns:
            bool: True if the compounds are supported, False otherwise.
        """

        num_supported = len([c for c in compounds if self.check_supported_compound(c) == True])

        if strict:
            # Checks all compounds are supported
            return num_supported == len(compounds)
        else:
            # Checks at least one is supported
            return num_supported > 0

    @abstractmethod
    def check_supported_compound(self, compound: str, strict: bool=True) -> bool:
        """
        Check if a single compound is supported by this package.
        Method can be overridden by subclasses to provide custom logic.
        Defalt behaviour assumes all "registered" compounds are supported.
        
        Args:
            compound (str): Name of the compound to check.
            strict (bool):  If True, all compounds must be supported.
                            If False, only a single compound must be supported.
        
        Returns:
            bool: True if the compound is supported, False otherwise.
        """
        return self.check_registered(compound)
    

class CompoundRegistry:

    def __init__(self):
        self.__compounds: Dict[str, Compound] = {}
        self.__packages: Dict[str, PropertyPackage] = {}
        self.__queue: Dict[str, list] = {"compounds": [], "packages": [], "bindings": []}
    
    @property
    def compounds(self):
        return self.__compounds
    
    @property
    def packages(self):
        return self.__packages
    
    def _build(self):
        # Building packages
        for package in self.__queue["packages"]:
            if isinstance(package, PropertyPackage):
                self._register_package(package)
            else:
                raise TypeError(f"Expected PropertyPackage, got {type(package)}")

        # Building compounds    
        for compound_name, sources, data in self.__queue["compounds"]:
            self._register_compound(compound_name, sources, data)

        # Building bindings
        for compound_name, package_name in self.__queue["bindings"]:
            if compound_name in self.compounds and package_name in self.packages:
                self._bind(compound_name, package_name)
            else:
                raise ValueError(f"Compound {compound_name} or Package {package_name} is not registered.")
        
        # Building dynamic bindings

    def queue_compound(self, compound_name: str, source: str, data: dict):
        """
        Queue a compound to be registered later.
        
        Args:
            compound_name (str): Name of the compound to queue.
            source (str): Name of the source providing the compound data.
            data (dict): Data associated with the compound from the source.
        """
        self.__queue["compounds"].append((compound_name, source, data))
    
    def queue_package(self, package: PropertyPackage):
        """
        Queue a property package to be registered later.
        
        Args:
            package (PropertyPackage): The property package to queue.
        """
        self.__queue["packages"].append(package)
    
    def queue_binding(self, compound_name: str, package_name: str):
        """
        Queue a binding between a compound and a property package.
        
        Args:
            compound_name (str): Name of the compound to bind.
            package_name (str): Name of the property package to bind to.
        """
        self.__queue["bindings"].append((compound_name, package_name))

    def _register_package(self, package: PropertyPackage):
        """
        Register a new property package in the registry.
        
        Args:
            package (PropertyPackage): The property package to register.
        """
        if package.name in self.packages:
            raise ValueError(f"Package {package.name} is already registered.")
        self.__packages[package.name] = package

    def _get_package(self, package_name: str) -> PropertyPackage:
        """
        Get a property package by its name.
        
        Args:
            package_name (str): Name of the property package to retrieve.
        
        Returns:
            PropertyPackage: The property package object if found, otherwise None.
        """
        return self.__packages.get(package_name, None)
    
    def _register_compound(self, compound: str, source: str, data: dict):
        """
        Register a new compound in the registry.
        
        Args:
            compound (Compound): Name of the compound to register.
            source (str): Name of the source providing the compound data.
            data (dict): Data associated with the compound from the source.
        """

        if not compound in self.compounds:
            # Create a new compound if it doesn't exist
            self.__compounds[compound] = Compound(compound)
        
        # Adding additional source to existing compound
        self.__compounds[compound].add_source(source, data)
    
    def _get_compound(self, compound_name: str) -> Compound:
        """
        Get a compound by its name.
        
        Args:
            compound_name (str): Name of the compound to retrieve.
        
        Returns:
            Compound: The compound object if found, otherwise None.
        """
        return self.__compounds.get(compound_name, None)
    
    def _bind(self, compound_name: str, package_name: str):
        """
        Bind a compound to a property package.

        Registers compound to package and vice versa.
        
        Args:
            compound_name (str): Name of the compound to bind.
            package_name (str): Name of the property package to bind to.
        
        Raises:
            ValueError: If either the compound or package is not registered.
        """
        if compound_name not in self.compounds:
            raise ValueError(f"Compound {compound_name} is not registered.")
        
        if package_name not in self.packages:
            raise ValueError(f"Package {package_name} is not registered.")
        
        # Add the package to the compound's list of supported packages
        self.__compounds[compound_name].add_package(package_name)

        # Add the compound to the package's list of supported compounds
        self.__packages[package_name].register_compound(compound_name)
    
    def _dynamic_bind(self, package_name: str):
        """
        Binds all compounds to 
        """
        if package_name not in self.packages:
            raise ValueError(f"Package {package_name} is not registered.")
        
        # Loop through all compounds
        # Check if compound is supported
        # Bind compound to package

    def _get_supported_packages(self, compounds: Set[str], strict=True) -> Set[str]:
        """
        Get a set of property packages that support the given compounds.
        
        Args:
            compounds (Set[str]): Set of compound names to check.
            strict (bool): If True, all compounds must be supported by all packages.
                           If False, at least one compound must be supported.
        
        Returns:
            Set[str]: Set of package names that support the given compounds.
        """

        supported_packages = set()
        
        # Looping through all registered packages
        for package in self.__packages.values():
            # Checking if package supports compounds with the given strictness
            if package.check_supported(compounds, strict):
                # If it does, add the package name to the supported packages
                supported_packages.add(package.name)
        
        return supported_packages
    
    def _get_supported_compounds(self, packages: Set[str], strict=True) -> Set[str]:
        """
        Get a set of compounds that are supported by the given property packages.
        
        Args:
            packages (Set[str]): Set of package names to check.
            strict (bool): If True, all packages must support all compounds.
                           If False, at least one package must support each compound.
        
        Returns:
            Set[str]: Set of compound names that are supported by the given packages.
        """

        supported_compounds = set()

        # Looping through all registered compounds
        
        # if strict:
        #     # Getting 
        #     for package in packages:
        #         p = self.get_package()
        #         if p is None:
        #             raise ValueError(f"Package {package} is not registered.")
        #         else:
        #             p.check_supported(compounds, strict)
        # else:
        # for package in packages:
        #     if 
        
        # return supported_compounds
    
    def _search_compounds(self,
            query: str, 
            package_filters: Dict[str]=None,
            filter_strict: bool=False
        ) -> set:

        """
        Searches for compounds based on the query.

        Args:
            query (str): The search query.
            package_filters (Dict[str], optional): Filters to apply on the compounds.
                If provided, only supported compounds from these packages will be considered.
            filter_strict (bool, optional): Only applies if package_filters is provided.
                                            If True, all packages must support all compounds.
                                            If False, at least one package must support each compound.
        
        Returns:
            list: List of compound names that match the query.
        
        """

        # Retrieving compounds by name
        filtered_compounds = Set(name for name in self.__compounds.keys() if query.lower() in name.lower())

        if package_filters is not None:
            filtered_compounds = self.get_supported_compounds(package_filters, strict=filter_strict).intersection(filtered_compounds)
        
        return filtered_compounds

    def _get_compound_names(self) -> list:
        """
        Returns a list of all compound names.

        Returns:
            list: List of compound names.
        """
        return list(self.__compounds.keys())


#  Simplified wrapper for back-end interactions
class RegistryView:
    def __init__(self, registry):
        self._registry = registry

    def register_package(self, pkg):
        self._registry.queue_package(pkg)

    def register_compound(self, name, source, data):
        self._registry.queue_compound(name, source, data)

    def bind(self, compound_name, package_name):
        self._registry.queue_bind(compound_name, package_name)


# Simplified wrapper for front-end integration
class RegistrySearch:

    def __init__(self, registry):
        self._registry = registry

    def search_compounds(self, query, package_filters=None, filter_strict=False):
        return self._registry.search_compounds(query, package_filters, filter_strict)

    def get_compound_names(self):
        return self._registry.get_compound_names()

    def get_supported_packages(self, compounds, strict=True):
        return self._registry.get_supported_packages(compounds, strict)
    
    def get_supported_compounds(self, packages, strict=True):
        return self._registry.get_supported_compounds(packages, strict)


registry = CompoundRegistry()
registry_view = RegistryView(registry)
registry_search = RegistrySearch(registry)

# Loading in loaders

from compounds.loaders import register_loader

# Building the registry

registry._build()

# export registry_search as registry <- need to do this

"""
"expected" workflow

registry = CompoundRegistry()

# Load in all seperate scripts registering to this registry

# Example workflow within a single loader script

    # Create PropertyPackage object with a unique name
    registry.register_package(PropertyPackage("name"))

    # Creates compound object in list, adds new source
    registry.register_compound("compound_name", "source_name", {"data": 123})

    # Default behaviour: linking compound to package
    registry.bind("compound_name", "package_name")

    # Alternative behaviour: generate links based on custom package
    registry.dynamic_bind("package_name")

registry._build()

"""

""" deprecated methods """

from compounds import deprecated

@deprecated
def get_compound(name: str) -> Compound | None:
    """
    Get a compound by its name.
    
    Args:
        name (str): Name of the compound to retrieve.
    
    Returns:
        compound (dict): the compound chemsep source data if found, otherwise None.
    """
    return registry._get_compound(name).get_source("chemsep")
