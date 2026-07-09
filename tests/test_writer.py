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
