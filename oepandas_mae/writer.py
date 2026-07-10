"""Maestro file writer for oepandas (``df.chem.to_mae``)."""
from __future__ import annotations

import contextlib
import logging
import os
import stat
import tempfile
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
from oepandas.arrays import MoleculeDtype
from openeye import oechem
from pandas.api.types import is_bool_dtype, is_float_dtype, is_integer_dtype, is_scalar

from ._compat import (
    WRITE_APPEND,
    WRITE_CREATE,
    OEMaestroTagConverter,
    OEMaestroWriter,
    OEMaestroWriterConfig,
    create_molecule_to_string_writer,
    is_gz,
)

log = logging.getLogger("oepandas-mae")

_DEFAULT_TITLE_COLUMN = "Title"

# One shared converter instance; IsFullMaestroKey is stateless.
_TAG_CONVERTER = OEMaestroTagConverter()


def _maestro_type_char(dtype) -> str:
    """Return the Maestro type character for a pandas dtype.

    :param dtype: A pandas/numpy dtype.
    :returns: ``"i"`` for integers, ``"r"`` for floats, ``"b"`` for booleans,
        ``"s"`` for everything else.
    """
    if is_bool_dtype(dtype):
        return "b"
    if is_integer_dtype(dtype):
        return "i"
    if is_float_dtype(dtype):
        return "r"
    return "s"


def _resolve_maestro_key(name: str, dtype, default_owner: str) -> tuple[str, str]:
    """Resolve a column name to a Maestro key and its effective type character.

    A name that is already a full Maestro key (``t_owner_name``) is used verbatim
    and its embedded type prefix is the effective type character. Otherwise the key
    is ``f"{type_char}_{default_owner}_{name}"`` where *type_char* is inferred from
    *dtype*.

    :param name: Column name.
    :param dtype: The column's pandas dtype (used only for constructed keys).
    :param default_owner: Owner segment for constructed keys.
    :returns: ``(maestro_key, effective_type_char)``.
    """
    if _TAG_CONVERTER.IsFullMaestroKey(name):
        return name, name[0]
    type_char = _maestro_type_char(dtype)
    return f"{type_char}_{default_owner}_{name}", type_char


def _format_value(value, type_char: str) -> str:
    """Format a cell value as Maestro SD-data text for a given type character.

    :param value: The cell value (assumed non-null; callers skip NaN/None).
    :param type_char: Effective type character (``"i"``/``"r"``/``"b"``/``"s"``).
    :returns: The textual representation to store as SD data.
    """
    if type_char == "b":
        if isinstance(value, str):
            token = value.strip().lower()
            if token in ("1", "true", "t", "yes", "y"):
                return "1"
            if token in ("0", "false", "f", "no", "n"):
                return "0"
            raise ValueError(f"Cannot interpret {value!r} as a boolean for a Maestro 'b_' property")
        return "1" if bool(value) else "0"
    if type_char == "i":
        return str(int(value))
    if type_char == "r":
        return str(float(value))
    return str(value)


def to_mae(
    self,
    fp: str | os.PathLike,
    molecule_column: str = "Molecule",
    *,
    title_column: str | None = "Title",
    columns: str | Iterable[str] | None = None,
    secondary_molecules_as: str | int = "smiles",
    secondary_molecule_flavor: int | None = None,
    index: bool = False,
    index_tag: str = "index",
    default_owner: str = "user",
    append: bool = False,
) -> None:
    """Write a DataFrame to a Maestro file. Registered as ``df.chem.to_mae``.

    Data columns are written as string SD data under computed Maestro keys, with the
    property type inferred from each column's dtype (see :func:`_resolve_maestro_key`).

    :param fp: Output path (``.mae``, ``.mae.gz``, ``.maegz``); gzip is auto-detected.
    :param molecule_column: Name of the primary molecule column (``MoleculeDtype``).
    :param title_column: Column used for the CT title, or ``None`` to skip. A missing
        default ``"Title"`` is skipped silently; a missing non-default name raises.
    :param columns: Data columns to write (``None`` = all non-molecule/non-title).
    :param secondary_molecules_as: Encoding for extra molecule columns (default SMILES).
    :param secondary_molecule_flavor: Optional flavor for secondary molecule encoding.
    :param index: If True, write the DataFrame index as a property.
    :param index_tag: Key used when ``index`` is True.
    :param default_owner: Owner segment for constructed (clean-name) keys.
    :param append: If True, append to an existing file (not allowed for gzip).
        **Note:** Append mode writes directly to the file and is NOT atomic: if writing
        fails partway through, previously written molecules from this call may remain
        in the file. Create mode (the default) is atomic — on failure the existing
        file is left unchanged.
    :returns: None.
    :raises KeyError: If the molecule column or a named column is absent.
    :raises TypeError: If the molecule column is not a ``MoleculeDtype``.
    :raises ValueError: On ``append`` with a gzip path, or a Maestro key collision.
    """
    df = self._obj

    # ---- Validate the primary molecule column ----
    if molecule_column not in df.columns:
        raise KeyError(f"Molecule column {molecule_column!r} not found in DataFrame")
    if not isinstance(df[molecule_column].dtype, MoleculeDtype):
        raise TypeError(f"Molecule column {molecule_column!r} is not a MoleculeDtype")

    # ---- Resolve title handling (value-based; see design Section 8) ----
    if (
        title_column is not None
        and title_column not in df.columns
        and title_column != _DEFAULT_TITLE_COLUMN
    ):
        raise KeyError(f"Title column {title_column!r} not found in DataFrame")
    write_title = title_column is not None and title_column in df.columns

    # ---- Resolve the columns to write ----
    if columns is None:
        excluded = {molecule_column}
        if title_column is not None:
            excluded.add(title_column)
        selected = [c for c in df.columns if c not in excluded]
    elif isinstance(columns, str):
        selected = [columns]
    else:
        selected = list(columns)
    for col in selected:
        if col not in df.columns:
            raise KeyError(f"Column {col!r} not found in DataFrame")

    # ---- Classify secondary molecule columns vs data columns ----
    secondary_columns = [
        c for c in selected
        if c != molecule_column and isinstance(df[c].dtype, MoleculeDtype)
    ]
    data_columns = [
        c for c in selected
        if c != molecule_column and c not in secondary_columns
    ]

    # ---- Resolve output keys and detect collisions before writing anything ----
    key_sources: dict[str, str] = {}

    def _register(key: str, source: str) -> None:
        if key in key_sources:
            raise ValueError(
                f"Maestro key collision: {source} and {key_sources[key]} "
                f"both map to {key!r}"
            )
        key_sources[key] = source

    data_keys: dict[str, tuple[str, str]] = {}
    for col in data_columns:
        key, type_char = _resolve_maestro_key(col, df[col].dtype, default_owner)
        _register(key, f"column {col!r}")
        data_keys[col] = (key, type_char)

    secondary_keys: dict[str, str] = {}
    for col in secondary_columns:
        # Resolve like data columns so full-key names pass through and all keys are
        # collision-checked. A clean molecule column (MoleculeDtype -> "s") yields
        # s_{default_owner}_{col} for its SMILES string.
        key, _ = _resolve_maestro_key(col, df[col].dtype, default_owner)
        _register(key, f"secondary molecule column {col!r}")
        secondary_keys[col] = key

    index_key: str | None = None
    if index:
        index_key = (
            index_tag if _TAG_CONVERTER.IsFullMaestroKey(index_tag)
            else f"s_{default_owner}_{index_tag}"
        )
        _register(index_key, "index")

    # ---- Open the writer ----
    fp_str = str(fp)
    if append and is_gz(Path(fp_str)):
        raise ValueError("append=True is not supported for gzip-compressed Maestro files")

    # ---- Symlink compression-mismatch guard (applies to both append and create) ----
    real_target = os.path.realpath(fp_str)
    if is_gz(Path(fp_str)) != is_gz(Path(real_target)):
        raise ValueError(
            f"Output path {fp_str!r} resolves to {real_target!r} with a different "
            f"compression format; refusing to write to avoid corrupting the target. "
            f"Pass a path whose extension matches the symlink target."
        )

    config = OEMaestroWriterConfig()
    config.SetMode(WRITE_APPEND if append else WRITE_CREATE)

    secondary_flavor = (
        secondary_molecule_flavor if isinstance(secondary_molecule_flavor, int) else None
    )
    secondary_to_string = create_molecule_to_string_writer(
        fmt=secondary_molecules_as,
        flavor=secondary_flavor,
        gzip=False,
        b64encode=False,
        strip=True,
    )

    # ---- Write rows ----
    def _write_rows(dest_path):
        with OEMaestroWriter(dest_path, config) as writer:
            for row_index, row in df.iterrows():
                mol = row[molecule_column]
                if isinstance(mol, oechem.OEMol):
                    out = oechem.OEGraphMol(mol) if mol.NumConfs() == 1 else oechem.OEMol(mol)
                elif isinstance(mol, oechem.OEMolBase):
                    out = oechem.OEGraphMol(mol)
                else:
                    out = oechem.OEGraphMol()

                # Output properties are driven solely by the DataFrame; drop any SD data
                # already on the source molecule so omitted/NaN columns cannot leak.
                oechem.OEClearSDData(out)

                if write_title:
                    out.SetTitle(str(row[title_column]))

                if index_key is not None:
                    oechem.OESetSDData(out, index_key, str(row_index))

                for col in secondary_columns:
                    value = row[col]
                    if value is not None:
                        oechem.OESetSDData(out, secondary_keys[col], secondary_to_string(value))

                for col in data_columns:
                    value = row[col]
                    # Scalar-safe: pd.isna on a list-like returns an array, so guard first.
                    if is_scalar(value) and pd.isna(value):
                        continue
                    key, type_char = data_keys[col]
                    oechem.OESetSDData(out, key, _format_value(value, type_char))

                if not writer.write(out):
                    raise OSError(f"Failed to write molecule at row {row_index} to {fp_str!r}")

    if append:
        _write_rows(fp_str)
    else:
        # Atomic overwrite; real_target already computed for symlink guard above.
        if fp_str.endswith(".mae.gz"):
            suffix = ".mae.gz"
        elif fp_str.endswith(".maegz"):
            suffix = ".maegz"
        else:
            suffix = ".mae"
        target_dir = os.path.dirname(real_target) or "."
        fd, tmp_path = tempfile.mkstemp(dir=target_dir, prefix="." + os.path.basename(real_target) + ".", suffix=suffix)
        os.close(fd)
        try:
            _write_rows(tmp_path)
        except BaseException:
            with contextlib.suppress(OSError):
                os.remove(tmp_path)
            raise
        if os.path.exists(real_target):
            os.chmod(tmp_path, stat.S_IMODE(os.stat(real_target).st_mode))
        else:
            current_umask = os.umask(0o022)
            os.umask(current_umask)
            os.chmod(tmp_path, 0o666 & ~current_umask)
        os.replace(tmp_path, real_target)
