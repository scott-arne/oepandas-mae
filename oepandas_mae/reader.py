"""Maestro file reader for oepandas."""
from __future__ import annotations

import logging
import os
from collections.abc import Iterable
from typing import Literal

import pandas as pd
from oepandas.arrays import MoleculeArray, MoleculeDtype
from openeye import oechem

from ._compat import (
    TAG_NAME,
    TAG_NONE,
    TAG_OWNER,
    TAG_TYPE,
    Dataset,
    MaestroMol,
    MaestroReader,
    MolConverter,
    OEMaestroReaderConfig,
    _add_smiles_columns,
)

log = logging.getLogger("oepandas-mae")

# Conformer test name -> OEChem conf test factory
_CONF_TESTS = {
    "default": oechem.OEDefaultConfTest,
    "absolute": oechem.OEAbsoluteConfTest,
    "absolute_canonical": oechem.OEAbsCanonicalConfTest,
    "isomeric": oechem.OEIsomericConfTest,
    "omega": oechem.OEOmegaConfTest,
}


def _format_tag_name(key: str, tag_flags: int) -> str:
    """Format a Maestro property key according to the tag bitmask.

    Maestro keys follow the ``t_o_name`` convention where *t* is the type
    character, *o* is the owner, and *name* is the property name.  The
    *tag_flags* bitmask controls which parts are included in the output.

    :param key: Raw Maestro property key (e.g. ``r_pdb_PDB_CRYST1_a``).
    :param tag_flags: Bitmask of TAG_TYPE, TAG_OWNER, TAG_NAME.
    :returns: Formatted tag name.
    """
    if len(key) < 4 or key[1] != "_":
        return key

    type_char = key[0]
    second_underscore = key.find("_", 2)
    if second_underscore == -1:
        return key

    owner = key[2:second_underscore]
    name = key[second_underscore + 1:]

    parts: list[str] = []
    if tag_flags & TAG_TYPE:
        parts.append(type_char)
    if tag_flags & TAG_OWNER:
        parts.append(owner)
    if tag_flags & TAG_NAME:
        parts.append(name)
    return "_".join(parts) if parts else key


def _read_maestro_data(
    filepath: str, config: OEMaestroReaderConfig
) -> tuple[list[oechem.OEGraphMol], list[dict[str, str]]]:
    """Read molecules and CT properties from a Maestro file.

    Uses the low-level ``MaestroReader`` and ``MolConverter`` to work around
    the SWIG binding limitation where CT properties are not transferred as SD
    data on the resulting OEGraphMol.

    :param filepath: Path to a Maestro file.
    :param config: Reader configuration (tags and perception).
    :returns: Tuple of (list of OEGraphMol, list of formatted CT property dicts).
    """
    try:
        reader = MaestroReader(filepath)
    except (FileNotFoundError, PermissionError):
        raise
    except Exception:
        log.debug("Failed to create MaestroReader for %s, returning empty lists", filepath)
        return [], []

    converter = MolConverter()
    converter.SetPerception(config.perception)

    mols: list[oechem.OEGraphMol] = []
    ct_props_list: list[dict[str, str]] = []

    mmol = MaestroMol()
    while reader.Read(mmol):
        mol = oechem.OEGraphMol()
        converter.Convert(mmol, mol)
        mols.append(mol)

        # Extract CT properties with tag formatting applied in Python
        formatted_props: dict[str, str] = {}
        if config.tags != TAG_NONE:
            for key in mmol.ct_properties:
                formatted_key = _format_tag_name(key, config.tags)
                formatted_props[formatted_key] = mmol.ct_properties[key]
        ct_props_list.append(formatted_props)

        mmol = MaestroMol()

    return mols, ct_props_list


def _group_conformers(
    mols: list[oechem.OEGraphMol],
    ct_props: list[dict[str, str]],
    conf_test,
) -> tuple[list[oechem.OEMol], list[dict[str, str]]]:
    """Group consecutive single-conformer molecules into multi-conformer OEMol objects.

    Uses the conf test's CompareMols method to decide if consecutive molecules
    should be grouped, and OEMol.NewConf to add conformers.  CT properties from
    the first molecule in each group are kept.

    :param mols: List of OEGraphMol from MaestroReader.
    :param ct_props: Parallel list of CT property dicts for each molecule.
    :param conf_test: OEChem conformer test object (e.g., OEDefaultConfTest()).
    :returns: Tuple of (grouped OEMol list, grouped CT properties list).
    """
    result_mols: list[oechem.OEMol] = []
    result_props: list[dict[str, str]] = []
    current: oechem.OEMol | None = None
    current_props: dict[str, str] | None = None

    for graph_mol, props in zip(mols, ct_props, strict=True):
        mol = oechem.OEMol(graph_mol)

        if current is None:
            current = mol
            current_props = props
            continue

        if conf_test.CompareMols(current, mol):
            current.NewConf(mol)
        else:
            result_mols.append(current)
            result_props.append(current_props)
            current = mol
            current_props = props

    if current is not None:
        result_mols.append(current)
        result_props.append(current_props)

    return result_mols, result_props


def _resolve_config(
    tags: int | None,
    perception: int | None,
    config: OEMaestroReaderConfig | None,
) -> OEMaestroReaderConfig:
    """Resolve an OEMaestroReaderConfig from kwargs and optional config object.

    :param tags: OEMaestroTag bitmask override.
    :param perception: OEMaestroPerception bitmask override.
    :param config: Base config object.
    :returns: Resolved OEMaestroReaderConfig.
    """
    if config is not None:
        resolved = OEMaestroReaderConfig()
        resolved.tags = config.tags
        resolved.perception = config.perception
    else:
        resolved = OEMaestroReaderConfig()
        # Default to TAG_NAME for clean column names (C++ default is TAG_ALL)
        resolved.tags = TAG_NAME

    # Kwargs override config
    if tags is not None:
        resolved.tags = tags
    if perception is not None:
        resolved.perception = perception

    return resolved


def read_mae(
    filepath: str | os.PathLike,
    *,
    molecule_column: str = "Molecule",
    title_column: str | None = "Title",
    no_title: bool = False,
    add_smiles: None | bool | str | Iterable[str] = None,
    usecols: None | str | Iterable[str] = None,
    numeric: None | str | dict[str, Literal["integer", "signed", "unsigned", "float"] | None] | Iterable[str] = None,
    conformer_test: Literal["default", "absolute", "absolute_canonical", "isomeric", "omega"] = "default",
    tags: int | None = None,
    perception: int | None = None,
    config: OEMaestroReaderConfig | None = None,
) -> pd.DataFrame:
    """Read structures from a Maestro file into a DataFrame.

    Reads .mae, .mae.gz, and .maegz files using oemaestro's native parser and
    produces a DataFrame following oepandas conventions.

    Use conformer_test to combine single conformers into multi-conformer molecules:
        - "default":
                Default conformer testing (groups consecutive molecules with matching titles).
        - "absolute":
                Combine conformers if they (1) have the same number of atoms and bonds in the
                same order, (2) each atom and bond have identical properties in the connection
                table, (3) have the same title.
        - "absolute_canonical":
                Combine conformers if they have the same canonical SMILES.
        - "isomeric":
                Combine conformers if they (1) have the same number of atoms and bonds in the
                same order, (2) each atom and bond have identical properties in the connection
                table, (3) have the same atom and bond stereochemistry, (4) have the same title.
        - "omega":
                Equivalent to "isomeric" except that invertible nitrogen stereochemistry is also
                taken into account.

    Use numeric to cast data columns to numeric types. SD data from Maestro files is always
    stored as strings:
        1. Single column name (default numeric cast)
        2. List of column names (default numeric cast)
        3. Dictionary of column names and specific numeric types to downcast to

    :param filepath: Path to a Maestro file (.mae, .mae.gz, .maegz).
    :param molecule_column: Name of the molecule column in the DataFrame.
    :param title_column: Name of the title column, or None to suppress.
    :param no_title: If True, suppress the title column.
    :param add_smiles: Include SMILES column(s) in the DataFrame.
    :param usecols: Subset of data tag columns to include.
    :param numeric: Data column(s) to make numeric.
    :param conformer_test: Conformer grouping strategy.
    :param tags: OEMaestroTag bitmask for column name formatting.
    :param perception: OEMaestroPerception bitmask.
    :param config: OEMaestroReaderConfig base configuration.
    :returns: Pandas DataFrame with molecule, title, and data columns.
    :raises ValueError: If conformer_test is not a valid option.
    """
    # ---- Validate conformer_test ----
    if conformer_test not in _CONF_TESTS:
        raise ValueError(
            f"Invalid conformer_test '{conformer_test}'. " f"Valid options: {list(_CONF_TESTS.keys())}"
        )

    # ---- Resolve title ----
    if no_title:
        title_column = None

    # ---- Resolve config ----
    resolved_config = _resolve_config(tags, perception, config)

    # ---- Normalize usecols ----
    if usecols is not None:
        usecols = frozenset((usecols,)) if isinstance(usecols, str) else frozenset(usecols)

    # ---- Normalize numeric ----
    if numeric is not None:
        if isinstance(numeric, str):
            numeric = {numeric: None}
        elif isinstance(numeric, Iterable) and not isinstance(numeric, dict):
            numeric = {col: None for col in numeric}
        if molecule_column in numeric:
            raise KeyError(f"Cannot make molecule column {molecule_column} numeric")

    # ---- Read molecules and CT properties ----
    filepath_str = str(filepath)
    raw_mols, ct_props_list = _read_maestro_data(filepath_str, resolved_config)

    # ---- Conformer grouping ----
    conf_test_obj = _CONF_TESTS[conformer_test]()
    mols_list, props_list = _group_conformers(raw_mols, ct_props_list, conf_test_obj)

    # ---- Handle empty file ----
    if len(mols_list) == 0:
        data = {molecule_column: pd.Series(dtype=MoleculeDtype())}
        if title_column is not None:
            data[title_column] = pd.Series(dtype=str)
        return pd.DataFrame(data)

    # ---- Build MoleculeArray ----
    mols = MoleculeArray(mols_list)

    # ---- Assemble primary columns ----
    data: dict[str, pd.Series] = {molecule_column: pd.Series(data=mols, dtype=MoleculeDtype())}

    if title_column is not None:
        data[title_column] = pd.Series(
            [mol.GetTitle() if isinstance(mol, oechem.OEMolBase) else None for mol in mols],
            dtype=str,
        )

    # ---- Extract CT properties into Dataset ----
    dataset = Dataset(usecols=usecols)
    for idx, props in enumerate(props_list):
        for tag, value in props.items():
            dataset.add(tag, value, idx)

    data = {**data, **dataset.to_series_dict()}

    # ---- Create DataFrame ----
    df = pd.DataFrame(data)

    # ---- Post-processing (only if non-empty) ----
    if len(df) > 0:
        # Numeric conversion
        if numeric is not None:
            for col, dtype in numeric.items():
                if isinstance(dtype, float):
                    dtype = "float"
                elif isinstance(dtype, int):
                    dtype = "integer"
                if col in df.columns:
                    try:
                        df[col] = pd.to_numeric(df[col], errors="coerce", downcast=dtype)
                    except Exception as e:
                        log.debug("Numeric downcast failed for column %s: %s", col, e)
                else:
                    log.warning("Column not found in DataFrame: %s", col)

        # SMILES columns
        if add_smiles is not None:
            _add_smiles_columns(df, molecule_column, add_smiles)

    return df
