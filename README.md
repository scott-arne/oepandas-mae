# oepandas-mae

Maestro file format support for [oepandas](https://github.com/scott-arne/oepandas). Reads `.mae`, `.mae.gz`, and
`.maegz` files into pandas DataFrames following oepandas conventions.

**This only requires a valid OpenEye Toolkit license. You do not need to have any Schrodinger software installed.**

## Installation

Install as an oepandas extra:

```bash
pip install oepandas[mae]
```

Or install standalone:

```bash
pip install oepandas-mae
```

**Requirements:** Python 3.10+, oepandas >= 3.2.0, oemaestro >= 0.1.0, OpenEye Toolkits

## Quick Start

```python
import oepandas as oepd

df = oepd.read_mae("tests/assets/5.maegz")
print(df)
```

Prints:

```text
                                            Molecule          Title ...
0  <oechem.OEMol; proxy of <Swig Object of type '...        Aspirin ...
1  <oechem.OEMol; proxy of <Swig Object of type '...      Ibuprofen ...
2  <oechem.OEMol; proxy of <Swig Object of type '...  Acetaminophen ...
3  <oechem.OEMol; proxy of <Swig Object of type '...       Caffeine ...
4  <oechem.OEMol; proxy of <Swig Object of type '...       Diazepam ...

[5 rows x 16 columns]
```

## Usage

### Column Naming

Rename the molecule and title columns, or suppress the title column entirely:

```python
import oepandas as oepd
df = oepd.read_mae("tests/assets/5.maegz", molecule_column="Mol", title_column="Name")
df = oepd.read_mae("tests/assets/5.maegz", no_title=True)
```

### Selecting Columns

Use `usecols` to include only specific data columns:

```python
df = read_mae("tests/assets/5.maegz", usecols=["NumAcceptors", "NumDonors"])
```

### Numeric Conversion

CT property values are strings by default. Use `numeric` to convert columns:

```python
# Convert a single column
df = read_mae("file.mae", numeric="pdb_tfactor")

# Convert multiple columns
df = read_mae("file.mae", numeric=["pdb_tfactor", "occupancy"])

# Specify downcast types
df = read_mae("file.mae", numeric={"pdb_tfactor": "float", "atom_count": "integer"})
```

### SMILES Columns

Add SMILES string columns alongside the molecule objects:

```python
df = read_mae("file.mae", add_smiles=True)
# Creates a "Molecule SMILES" column
```

### Conformer Grouping

Consecutive structures in a Maestro file can be grouped into multi-conformer molecules:

```python
df = read_mae("file.mae", conformer_test="absolute")
```

Available conformer tests:
- `"default"` -- groups consecutive molecules with matching titles
- `"absolute"` -- requires identical atom/bond ordering, properties, and title
- `"absolute_canonical"` -- requires matching canonical SMILES
- `"isomeric"` -- like absolute, but also requires matching stereochemistry
- `"omega"` -- like isomeric, plus invertible nitrogen stereochemistry

### Tag Formatting

Maestro property keys follow a `type_owner_name` convention (e.g., `r_pdb_PDB_CRYST1_a`). By default, only the name portion is used as the column name. Control this with `tags`:

```python
from oepandas_mae import read_mae, TAG_ALL, TAG_NAME, TAG_NONE

# Default: clean names only (e.g., "PDB_CRYST1_a")
df = read_mae("file.mae")

# Full Maestro keys (e.g., "r_pdb_PDB_CRYST1_a")
df = read_mae("file.mae", tags=TAG_ALL)

# No data columns (molecules and titles only)
df = read_mae("file.mae", tags=TAG_NONE)
```

### Perception Control

Control post-parse chemical perception with the `perception` parameter:

```python
from oepandas_mae import read_mae, PERCEPTION_NONE

df = read_mae("file.mae", perception=PERCEPTION_NONE)
```

### Configuration Object

For repeated use, pass an `OEMaestroReaderConfig` object. Keyword arguments override config values:

```python
from oepandas_mae import read_mae, OEMaestroReaderConfig, TAG_ALL

config = OEMaestroReaderConfig()
config.tags = TAG_ALL

df = read_mae("file.mae", config=config)
df = read_mae("file.mae", config=config, tags=TAG_NAME)  # tags overrides config
```

## License

MIT
