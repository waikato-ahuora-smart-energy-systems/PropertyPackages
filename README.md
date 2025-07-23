To install this:

```sh
pip install git+https://github.com/waikato-ahuora-smart-energy-systems/PropertyPackages.git
```

To install the development version, after cloning this repo:

```sh
pip install -e .
```

To use

# About

This repository is split into two main tools: the compound database and property package builder.

# Compound Database

## Usage

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

### Depreciated methods

#### > get_compound(name)

> [!WARNING]
> method is depreciated

```python
# example usage
from compounds import CompoundDB as db
benzene = db.get_compound("benzene")
benzene["AbsEntropy"].value # "269300"
benzene["AbsEntropy"].unit # "J/kmol/K" 
```

#### > get_compound_names

> [!WARNING]
> method is depreciated

```python
# example usage
from compounds import CompoundDB as db
names = db.get_compound_names()
print(names)
> ["benzene", "methane", "butane", "..."]
```

#### > search_compounds(query)

> [!WARNING]
> method is depreciated

```python
# example usage
from compounds import CompoundDB as db
names = db.search_compounds("hex")
print(names)
> ["hexane", "1-hexane", "..."]
```

## Next Steps



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