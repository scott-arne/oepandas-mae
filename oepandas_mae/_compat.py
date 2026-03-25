"""Compatibility layer for oepandas internal imports and oemaestro re-exports.

Centralizes dependencies on oepandas private API so there is exactly one place
to update if oepandas refactors its internals.
"""
from oepandas.pandas_extensions import Dataset, _add_smiles_columns

from oemaestro import (
    TAG_NONE,
    TAG_TYPE,
    TAG_OWNER,
    TAG_NAME,
    TAG_ALL,
    PERCEPTION_NONE,
    PERCEPTION_CONNECTIVITY,
    PERCEPTION_RINGS,
    PERCEPTION_BOND_ORDERS,
    PERCEPTION_IMPLICIT_HYDROGENS,
    PERCEPTION_FORMAL_CHARGES,
    PERCEPTION_ALL,
    OEMaestroReaderConfig,
    OEMaestroReader,
)

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
    "PERCEPTION_ALL",
    "OEMaestroReaderConfig",
    "OEMaestroReader",
]
