import pandas as pd


def test_writer_compat_symbols_importable():
    from oepandas_mae._compat import (
        WRITE_APPEND,
        WRITE_CREATE,
        OEMaestroTagConverter,
        OEMaestroWriter,
        OEMaestroWriterConfig,
        create_molecule_to_string_writer,
        is_gz,
    )

    assert OEMaestroWriter is not None
    assert OEMaestroWriterConfig is not None
    assert OEMaestroTagConverter is not None
    assert WRITE_CREATE != WRITE_APPEND
    assert callable(create_molecule_to_string_writer)
    assert callable(is_gz)


def test_maestro_type_char():
    from oepandas_mae.writer import _maestro_type_char

    assert _maestro_type_char(pd.Series([1], dtype="int64").dtype) == "i"
    assert _maestro_type_char(pd.Series([1.0], dtype="float64").dtype) == "r"
    assert _maestro_type_char(pd.Series([True], dtype="bool").dtype) == "b"
    assert _maestro_type_char(pd.Series(["x"], dtype="object").dtype) == "s"


def test_resolve_maestro_key_clean_name_infers_type():
    from oepandas_mae.writer import _resolve_maestro_key

    key, type_char = _resolve_maestro_key("count", pd.Series([1], dtype="int64").dtype, "user")
    assert key == "i_user_count"
    assert type_char == "i"


def test_resolve_maestro_key_full_key_passthrough():
    from oepandas_mae.writer import _resolve_maestro_key

    # Embedded prefix wins over dtype inference; object dtype would otherwise be "s".
    key, type_char = _resolve_maestro_key(
        "r_m_pdb_tfactor", pd.Series(["x"], dtype="object").dtype, "user"
    )
    assert key == "r_m_pdb_tfactor"
    assert type_char == "r"


def test_resolve_maestro_key_respects_default_owner():
    from oepandas_mae.writer import _resolve_maestro_key

    key, _ = _resolve_maestro_key("score", pd.Series([1.0], dtype="float64").dtype, "myapp")
    assert key == "r_myapp_score"


def test_format_value():
    from oepandas_mae.writer import _format_value

    assert _format_value(42, "i") == "42"
    assert _format_value(3.5, "r") == "3.5"
    assert _format_value(True, "b") == "1"
    assert _format_value(False, "b") == "0"
    assert _format_value("hello", "s") == "hello"
