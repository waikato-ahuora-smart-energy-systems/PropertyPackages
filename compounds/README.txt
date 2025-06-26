# Ahuora Compounds Database

## Loading and using data

### Simple example workflow

```python

from compounds import CompoundDB as db

benzene = db.get_compound("benzene")

benzene_cp = benzene.get_source("chemsep")

benzene_cp["AbsEntropy"].value # "269300"
benzene_cp["AbsEntropy"].unit # "J/kmol/K"

```

### Methods

#### Active



#### Depreciated

get_compound

get_compound_names

get_compounds

## Adding new information

