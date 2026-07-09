"""Compatibility layer for oepandas internal imports and oemaestro re-exports.

Centralizes dependencies on oepandas private API so there is exactly one place
to update if oepandas refactors its internals.
"""
from oemaestro import (
    PERCEPTION_ALL,
    PERCEPTION_BOND_ORDERS,
    PERCEPTION_CONNECTIVITY,
    PERCEPTION_DEFAULT,
    PERCEPTION_FORMAL_CHARGES,
    PERCEPTION_IMPLICIT_HYDROGENS,
    PERCEPTION_NONE,
    PERCEPTION_RINGS,
    TAG_ALL,
    TAG_NAME,
    TAG_NONE,
    TAG_OWNER,
    TAG_TYPE,
    WRITE_APPEND,
    WRITE_CREATE,
    MaestroMol,
    MaestroReader,
    MolConverter,
    OEMaestroReader,
    OEMaestroReaderConfig,
    OEMaestroTagConverter,
    OEMaestroWriter,
    OEMaestroWriterConfig,
)
from oepandas.pandas_extensions import Dataset, _add_smiles_columns
from oepandas.util import create_molecule_to_string_writer, is_gz

__all__ = [
    "Dataset",
    "_add_smiles_columns",
    "TAG_NONE",
    "TAG_TYPE",
    "TAG_OWNER",
    "TAG_NAME",
    "TAG_ALL",
    "PERCEPTION_NONE",
    "PERCEPTION_CONNECTIVITY",
    "PERCEPTION_RINGS",
    "PERCEPTION_BOND_ORDERS",
    "PERCEPTION_IMPLICIT_HYDROGENS",
    "PERCEPTION_FORMAL_CHARGES",
    "PERCEPTION_DEFAULT",
    "PERCEPTION_ALL",
    "OEMaestroReaderConfig",
    "OEMaestroReader",
    "MaestroReader",
    "MaestroMol",
    "MolConverter",
    "WRITE_CREATE",
    "WRITE_APPEND",
    "OEMaestroWriter",
    "OEMaestroWriterConfig",
    "OEMaestroTagConverter",
    "create_molecule_to_string_writer",
    "is_gz",
]
