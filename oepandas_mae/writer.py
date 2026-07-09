"""Maestro file writer for oepandas (``df.chem.to_mae``)."""
from __future__ import annotations

from pandas.api.types import is_bool_dtype, is_float_dtype, is_integer_dtype

from ._compat import OEMaestroTagConverter

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
        return "1" if bool(value) else "0"
    if type_char == "i":
        return str(int(value))
    if type_char == "r":
        return str(float(value))
    return str(value)
