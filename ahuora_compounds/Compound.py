from typing import Dict, Set

class Compound:
    def __init__(self, name: str):
        self.__name = name
        self.__sources: Dict[str, dict] = {}
        self.__packages: Set[str] = set()  # supported packages

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
