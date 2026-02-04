from typing import Dict, Set
import ahuora_compounds.loaders as loaders
import pkgutil
import importlib
from ahuora_compounds.Compound import Compound
from ahuora_compounds.PropertyPackage import PropertyPackage
from ahuora_compounds.RegistryLoader import RegistryLoader


class CompoundRegistry:

    def __init__(self):
        self.__compounds: Dict[str, Compound] = {}
        self.__packages: Dict[str, PropertyPackage] = {}
        self.__queue: Dict[str, list] = {"compounds": [], "packages": [], "bindings": [], "dynamic_bindings": []}
        self._built: bool = False

    @property
    def compounds(self):
        return self.__compounds

    @property
    def packages(self):
        return self.__packages

    def _build(self):
        if not self._built:
            from ahuora_compounds.loaders import loaders_list

            self._built = True

            # Import all loaders
            for loader in loaders_list:
                loader(RegistryLoader(self))

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
            for package_name in self.__queue["dynamic_bindings"]:
                if package_name in self.packages:
                    self._dynamic_bind(package_name)
                else:
                    raise ValueError(f"Package {package_name} is not registered.")
    
    def _discover_loaders(self):
        # Loading modules
        for module in pkgutil.iter_modules(loaders.__path__):
            #  class pkgutil.ModuleInfo(module_finder, name, ispkg)
            importlib.import_module(f"compounds.loaders.{module.name}")

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
    
    def queue_dynamic_binding(self, package_name: str):
        """
        Queue a dynamic binding for a property package.
        
        Args:
            package_name (str): Name of the property package to bind dynamically.
        """
        self.__queue["dynamic_bindings"].append(package_name)

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

        if compound not in self.compounds:

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
        Binds all compounds to a package dynamically. (based on package check_supported_compounds implementation)
        """
        if package_name not in self.packages:
            raise ValueError(f"Package {package_name} is not registered.")

        # Loop through all compounds
        for compound_name in self.compounds:
            # Check if compound is supported by package
            if self.__packages[package_name].check_supported_compound(compound_name):
                # Bind compound to package
                self._bind(compound_name, package_name)

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
            if package.check_supported_compounds(compounds, strict):
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
        for package in packages:
            p = self._get_package(package)
            if p is None:
                raise ValueError(f"Package {package} is not registered.")
            else:
                if strict:
                    if len(supported_compounds) == 0:
                        supported_compounds = p.compounds
                    else:
                        supported_compounds = supported_compounds.intersection(p.compounds)     
                else:
                    supported_compounds = supported_compounds.union(p.compounds)
        
        return supported_compounds

    def _search_compounds(self,
                          query: str,
                          package_filters: Set[str] = None,
                          filter_strict: bool = False
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
        filtered_compounds = [name for name in self.__compounds.keys() if query.lower() in name.lower()]

        if package_filters is not None:
            filtered_compounds = self._get_supported_compounds(package_filters, strict=filter_strict) \
                .intersection(filtered_compounds)

        return filtered_compounds

    def _get_compound_names(self) -> list:
        """
        Returns a list of all compound names.

        Returns:
            list: List of compound names.
        """
        return list(self.__compounds.keys())

    # To force build module on import
    # https://stackoverflow.com/questions/77186124/run-a-function-every-time-a-method-in-a-class-is-called
    def __getattribute__(self, attr):
        # Retrieve the attribute
        attr = super().__getattribute__(attr)
        # Check if attribute is a method
        if callable(attr):
            # Ensure build is after discovery
            if not attr.__name__ == '_discover_loaders':
                super().__getattribute__('_build')()
        return attr

        
