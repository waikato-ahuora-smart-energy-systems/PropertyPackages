## About

This repository is split into two main tools: the compound database and property package builder.

## Installation

To install this:

```sh
pip install git+https://github.com/waikato-ahuora-smart-energy-systems/PropertyPackages.git
```

To install the development version, after cloning this repo:

```sh
pip install -e .
```

To use



## Compound Database

### Usage

```python
# example usage
from compounds.CompoundDB import db

db.search_compounds("ane")
> ["octane", "hexane", "..."]

db.get_compound_names()
> ["amonia", "benzene", "carbon dioxide", "..."]

db.get_compound(name)
> <Compound Object>

db.get_supported_packages(["benzene", "toluene", "..."])
> ["peng-robinson", "helmholtz"]

db.get_supported_compounds(["peng-robinson", "helmholtz"])
> ["benzene", "toluene", "..."]
```

### Development

```

Creating new data loader under compounds/loaders/

3 step process

1 - register package
2 - register compounds
3 - bind compounds to package

```

```python

# general workflow

# default property package

@loader("example_loader")
def load(registry):

    registry.register_package(DefaultPropertyPackage("example_pp"))
    
    registry.register_compound('compound1', "example_pp", {})
    registry.register_compound('compound2', "example_pp", {})

    registry.bind("compound1", "example_pp")
    registry.bind("compound2", "example_pp")

# custom property package

class MilkPropertyPackage(PropertyPackage):
    def check_supported_compound(self, compound: str, strict: bool = True) -> bool:
        return True

@loader("example_loader")
def load(registry):

    registry.register_package(MilkPropertyPackage("example_custom_pp"))
    
    registry.register_compound('compound1', "example_custom_pp", {})
    registry.register_compound('compound2', "example_custom_pp", {})
    registry.register_compound('compound3', "example_custom_pp", {})

    # Generates based on PropertyPackage check_supported_compound implementation
    # In this case all compounds across loaders will be valid for this property package
    registry.dynamic_bind("example_custom_pp") 

```

```
more examples in compounds/loaders/
```

### Deprecated methods

#### > get_compound(name)

> [!WARNING]
> method is deprecated

```python
# example usage
from CompoundDB import get_compound
benzene = get_compound("benzene")
benzene["AbsEntropy"].value # "269300"
benzene["AbsEntropy"].unit # "J/kmol/K" 
```

#### > get_compound_names

> [!WARNING]
> method is deprecated

```python
# example usage
from CompoundDB import get_compound_names
names = get_compound_names()
print(names)
> ["benzene", "methane", "butane", "..."]
```

#### > search_compounds(query)

> [!WARNING]
> method is deprecated

```python
# example usage
from compounds import search_compounds
names = search_compounds("hex")
print(names)
> ["hexane", "1-hexane", "..."]
```

# Property Package Builder

## Usage

```python
import property_packages
property_packages.build_package("helmholtz",["h2o"])
# Returns IDAES compatible ParameterBlock
```

## Next Steps

Implement something to translate between apis, e.g for if the user specifies in Temperature and mass flow but the property package uses flow_mol and enthalpy

```py
property_package.build_state(
    "helmholtz",
    "h2o",
    {
        "temperature",
        "pressure",
        "mass_flow",
    }
) -> {
    "enthalpy":
    "pressure":
    "quality":
    "flow_mol":
}
```

Support different phases:

```py
property_packages.build_package("helmholtz",["h2o"],["vap","liq"])
```