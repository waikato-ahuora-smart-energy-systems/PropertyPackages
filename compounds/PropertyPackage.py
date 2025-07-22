from abc import ABC, abstractmethod
from typing import Set

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

    def check_supported_compounds(self, compounds: Set[str], strict: bool = True):
        """
        Check if a set of compounds are supported by this package.
        
        Args:
            compounds (Set[str]): Name of the compound to check.
            strict (bool):  If True, all compounds must be supported.
                            If False, only a single compound must be supported.
        
        Returns:
            bool: True if the compounds are supported, False otherwise.
        """

        num_supported = len([c for c in compounds if self.check_supported_compound(c)])

        if strict:
            # Checks all compounds are supported
            return num_supported == len(compounds)
        else:
            # Checks at least one is supported
            return num_supported > 0

    @abstractmethod
    def check_supported_compound(self, compound: str, strict: bool = True) -> bool:
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