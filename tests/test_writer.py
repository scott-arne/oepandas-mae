import pandas as pd
import pytest
from oepandas.arrays import MoleculeArray, MoleculeDtype
from openeye import oechem

from oepandas_mae import TAG_ALL
from oepandas_mae.reader import read_mae


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


def _make_mol(title="mol"):
    """Build a minimal 2-atom OEMol with coordinates and a title."""
    mol = oechem.OEMol()
    c = mol.NewAtom(6)
    o = mol.NewAtom(8)
    mol.NewBond(c, o, 2)
    mol.SetCoords(oechem.OEFloatArray([0.0, 0.0, 0.0, 1.2, 0.0, 0.0]))
    mol.SetTitle(title)
    return mol


def _mol_frame(mols, **columns):
    """Build a DataFrame with a MoleculeDtype 'Molecule' column plus extra columns."""
    data = {"Molecule": pd.Series(MoleculeArray(mols), dtype=MoleculeDtype())}
    data.update(columns)
    return pd.DataFrame(data)


def test_to_mae_registered_on_accessor():
    import oepandas_mae  # noqa: F401

    df = _mol_frame([_make_mol("a")])
    assert hasattr(df.chem, "to_mae")
    assert callable(df.chem.to_mae)


def test_write_read_roundtrip_structure(tmp_path):
    df = _mol_frame([_make_mol("aspirin"), _make_mol("caffeine")])
    out = tmp_path / "out.mae"
    df.chem.to_mae(out)

    back = read_mae(out)
    assert len(back) == 2
    assert list(back["Title"]) == ["aspirin", "caffeine"]
    assert all(m.NumAtoms() == 2 for m in back["Molecule"])


def test_missing_molecule_column_raises(tmp_path):
    df = pd.DataFrame({"NotMol": [1, 2]})
    with pytest.raises(KeyError, match="not found"):
        df.chem.to_mae(tmp_path / "x.mae")


def test_non_molecule_column_raises(tmp_path):
    df = pd.DataFrame({"Molecule": [1, 2]})
    with pytest.raises(TypeError, match="MoleculeDtype"):
        df.chem.to_mae(tmp_path / "x.mae")


def test_type_preservation_roundtrip(tmp_path):
    df = _mol_frame(
        [_make_mol("m")],
        count=pd.array([42], dtype="int64"),
        energy=pd.array([3.5], dtype="float64"),
        active=pd.array([True], dtype="bool"),
        name=pd.array(["hello"], dtype="object"),
    )
    out = tmp_path / "typed.mae"
    df.chem.to_mae(out)

    back = read_mae(out, tags=TAG_ALL)
    row = back.iloc[0]
    assert row["i_user_count"] == 42
    assert row["r_user_energy"] == 3.5
    assert int(row["b_user_active"]) == 1
    assert row["s_user_name"] == "hello"


def test_string_column_roundtrip(tmp_path):
    df = _mol_frame([_make_mol("m")], label=pd.array(["abc"], dtype="object"))
    out = tmp_path / "s.mae"
    df.chem.to_mae(out)

    back = read_mae(out)
    assert back.iloc[0]["label"] == "abc"


def test_columns_selection_filters_data_columns(tmp_path):
    df = _mol_frame(
        [_make_mol("m")],
        keep=pd.array([1], dtype="int64"),
        drop=pd.array([2], dtype="int64"),
    )
    out = tmp_path / "sel.mae"
    df.chem.to_mae(out, columns=["keep"])

    back = read_mae(out, tags=TAG_ALL)
    assert "i_user_keep" in back.columns
    assert "i_user_drop" not in back.columns


def test_default_owner_used_in_keys(tmp_path):
    df = _mol_frame([_make_mol("m")], score=pd.array([1.0], dtype="float64"))
    out = tmp_path / "own.mae"
    df.chem.to_mae(out, default_owner="myapp")

    back = read_mae(out, tags=TAG_ALL)
    assert "r_myapp_score" in back.columns


def test_title_column_used_when_present(tmp_path):
    df = _mol_frame([_make_mol("intrinsic")], Title=["override"])
    out = tmp_path / "t.mae"
    df.chem.to_mae(out)

    back = read_mae(out)
    assert back.iloc[0]["Title"] == "override"


def test_missing_default_title_uses_intrinsic(tmp_path):
    df = _mol_frame([_make_mol("intrinsic")])  # no "Title" column
    out = tmp_path / "t2.mae"
    df.chem.to_mae(out)  # default title_column="Title" absent -> skip silently

    back = read_mae(out)
    assert back.iloc[0]["Title"] == "intrinsic"


def test_missing_nondefault_title_raises(tmp_path):
    df = _mol_frame([_make_mol("m")])
    with pytest.raises(KeyError, match="Title column"):
        df.chem.to_mae(tmp_path / "t3.mae", title_column="Name")


def test_key_collision_raises(tmp_path):
    df = _mol_frame(
        [_make_mol("m")],
        count=pd.array([1], dtype="int64"),          # -> i_user_count
        **{"i_user_count": pd.array([2], dtype="int64")},  # full key -> i_user_count
    )
    with pytest.raises(ValueError, match="collision"):
        df.chem.to_mae(tmp_path / "c.mae")


def test_none_molecule_does_not_raise(tmp_path):
    df = _mol_frame([_make_mol("a"), None])
    out = tmp_path / "n.mae"
    df.chem.to_mae(out)  # None -> empty OEGraphMol; must not raise

    back = read_mae(out)
    assert "a" in list(back["Title"])


def test_empty_dataframe_writes_empty_file(tmp_path):
    empty = _mol_frame([])
    out = tmp_path / "empty.mae"
    empty.chem.to_mae(out)

    back = read_mae(out)
    assert len(back) == 0
    assert "Molecule" in back.columns


def test_append_mode_appends(tmp_path):
    out = tmp_path / "app.mae"
    _mol_frame([_make_mol("first")]).chem.to_mae(out)
    _mol_frame([_make_mol("second")]).chem.to_mae(out, append=True)

    back = read_mae(out)
    assert list(back["Title"]) == ["first", "second"]


def test_append_gzip_raises(tmp_path):
    df = _mol_frame([_make_mol("m")])
    with pytest.raises(ValueError, match="append"):
        df.chem.to_mae(tmp_path / "x.mae.gz", append=True)


def test_gzip_output_roundtrip(tmp_path):
    df = _mol_frame([_make_mol("g")], val=pd.array([7], dtype="int64"))
    out = tmp_path / "out.mae.gz"
    df.chem.to_mae(out)

    back = read_mae(out, tags=TAG_ALL)
    assert back.iloc[0]["i_user_val"] == 7


def test_multiconformer_roundtrip(tmp_path, conformers_mae):
    df = read_mae(conformers_mae, conformer_test="absolute")  # 1 grouped multi-conf row
    assert df.iloc[0]["Molecule"].NumConfs() == 2
    out = tmp_path / "conf.mae"
    df.chem.to_mae(out)

    back = read_mae(out, conformer_test="absolute")
    assert len(back) == 1
    assert back.iloc[0]["Molecule"].NumConfs() == 2


def test_fixture_roundtrips(tmp_path, simple_mae, multi_mae, protein_mae, protein_maegz, conformers_mae):
    for src in (simple_mae, multi_mae, protein_mae, protein_maegz, conformers_mae):
        df = read_mae(src, tags=TAG_ALL)  # full keys -> exact passthrough, no owner drift
        out = tmp_path / f"{src.stem}_out.mae"
        df.chem.to_mae(out)
        back = read_mae(out, tags=TAG_ALL)
        assert len(back) == len(df)
        assert list(back["Title"]) == list(df["Title"])
        assert [m.NumAtoms() for m in back["Molecule"]] == [m.NumAtoms() for m in df["Molecule"]]
        assert [m.NumBonds() for m in back["Molecule"]] == [m.NumBonds() for m in df["Molecule"]]


def test_coordinates_roundtrip(tmp_path, simple_mae):
    df = read_mae(simple_mae)
    mol = df.iloc[0]["Molecule"]
    original = {a.GetIdx(): tuple(round(c, 3) for c in mol.GetCoords(a)) for a in mol.GetAtoms()}

    out = tmp_path / "coords.mae"
    df.chem.to_mae(out)
    back = read_mae(out)
    back_mol = back.iloc[0]["Molecule"]
    for atom in back_mol.GetAtoms():
        x, y, z = back_mol.GetCoords(atom)
        ox, oy, oz = original[atom.GetIdx()]
        assert abs(x - ox) < 0.01
        assert abs(y - oy) < 0.01
        assert abs(z - oz) < 0.01


def test_custom_molecule_and_title_columns(tmp_path):
    df = _mol_frame([_make_mol("keep-title")])
    df = df.rename(columns={"Molecule": "Mol"})
    df["Name"] = ["renamed"]
    out = tmp_path / "custom.mae"
    df.chem.to_mae(out, molecule_column="Mol", title_column="Name")

    back = read_mae(out)
    assert back.iloc[0]["Title"] == "renamed"
    assert back.iloc[0]["Molecule"].NumAtoms() == 2


def test_columns_selection_single_string(tmp_path):
    df = _mol_frame(
        [_make_mol("m")],
        keep=pd.array([1], dtype="int64"),
        drop=pd.array([2], dtype="int64"),
    )
    out = tmp_path / "sel_str.mae"
    df.chem.to_mae(out, columns="keep")

    back = read_mae(out, tags=TAG_ALL)
    assert "i_user_keep" in back.columns
    assert "i_user_drop" not in back.columns


def test_object_listlike_value_does_not_raise(tmp_path):
    df = _mol_frame([_make_mol("m")], payload=pd.Series([[1, 2]], dtype="object"))
    out = tmp_path / "obj.mae"
    df.chem.to_mae(out)  # list-like object cell -> str(value); must not raise

    back = read_mae(out)
    assert back.iloc[0]["payload"] == "[1, 2]"


def test_preexisting_sd_data_does_not_leak(tmp_path):
    mol = _make_mol("m")
    oechem.OESetSDData(mol, "s_user_stale", "leaked")  # not from any DataFrame column
    df = _mol_frame([mol], keep=pd.array([1], dtype="int64"))
    out = tmp_path / "leak.mae"
    df.chem.to_mae(out)

    back = read_mae(out, tags=TAG_ALL)
    assert "i_user_keep" in back.columns
    assert "s_user_stale" not in back.columns  # stale SD data was cleared


def test_full_key_boolean_string_false_not_flipped(tmp_path):
    df = _mol_frame([_make_mol("m")], **{"b_user_flag": pd.array(["0"], dtype="object")})
    out = tmp_path / "boolstr.mae"
    df.chem.to_mae(out)
    back = read_mae(out, tags=TAG_ALL)
    # b_ reads back as int 0/1; "0" must stay 0, not flip to 1
    assert int(back.iloc[0]["b_user_flag"]) == 0


def test_full_key_boolean_ambiguous_string_raises(tmp_path):
    df = _mol_frame([_make_mol("m")], **{"b_user_flag": pd.array(["maybe"], dtype="object")})
    with pytest.raises(ValueError, match="boolean"):
        df.chem.to_mae(tmp_path / "boolbad.mae")


def test_write_failure_raises(tmp_path, monkeypatch):
    import oepandas_mae.writer as writer_mod

    class _FailWriter:
        def __init__(self, *a, **k): pass
        def write(self, mol): return False
        def close(self): pass
        def __enter__(self): return self
        def __exit__(self, *a): return False

    monkeypatch.setattr(writer_mod, "OEMaestroWriter", _FailWriter)
    df = _mol_frame([_make_mol("m")])
    with pytest.raises(OSError, match="Failed to write"):
        df.chem.to_mae(tmp_path / "fail.mae")


def test_create_mode_failure_preserves_existing_file(tmp_path):
    out = tmp_path / "existing.mae"
    _mol_frame([_make_mol("original")]).chem.to_mae(out)
    original_bytes = out.read_bytes()
    # An overwrite that fails mid-write (ambiguous boolean string) must not clobber the original.
    bad = _mol_frame([_make_mol("new")], **{"b_user_flag": pd.array(["maybe"], dtype="object")})
    with pytest.raises(ValueError):
        bad.chem.to_mae(out)
    assert out.read_bytes() == original_bytes
    # And no leftover temp files in the directory.
    assert [p.name for p in tmp_path.iterdir()] == ["existing.mae"]


def test_overwrite_preserves_file_permissions(tmp_path):
    import os
    import stat

    out = tmp_path / "shared.mae"
    _mol_frame([_make_mol("first")]).chem.to_mae(out)
    os.chmod(out, 0o644)
    _mol_frame([_make_mol("second")]).chem.to_mae(out)  # atomic overwrite
    mode = stat.S_IMODE(out.stat().st_mode)
    assert mode == 0o644, f"expected 0o644, got {oct(mode)}"
    # sanity: the overwrite actually happened
    from oepandas_mae.reader import read_mae
    assert list(read_mae(out)["Title"]) == ["second"]


def test_secondary_molecule_written_as_smiles(tmp_path):
    df = _mol_frame([_make_mol("primary")])
    df["Other"] = pd.Series(MoleculeArray([_make_mol("second")]), dtype=MoleculeDtype())
    out = tmp_path / "sec.mae"
    df.chem.to_mae(out)

    back = read_mae(out)
    assert "Other" in back.columns
    value = back.iloc[0]["Other"]
    assert isinstance(value, str) and value  # non-empty SMILES string


def test_index_written_as_property(tmp_path):
    df = _mol_frame([_make_mol("a"), _make_mol("b")])
    df.index = [10, 20]
    out = tmp_path / "idx.mae"
    df.chem.to_mae(out, index=True)

    back = read_mae(out)
    assert list(back["index"].astype(str)) == ["10", "20"]


def test_index_custom_tag(tmp_path):
    df = _mol_frame([_make_mol("a")])
    out = tmp_path / "idx2.mae"
    df.chem.to_mae(out, index=True, index_tag="row_id")

    back = read_mae(out)
    assert "row_id" in back.columns


def test_index_collision_with_column_raises(tmp_path):
    df = _mol_frame([_make_mol("m")], index=pd.array(["dup"], dtype="object"))
    # data column "index" -> s_user_index; index_tag "index" -> s_user_index -> collision
    with pytest.raises(ValueError, match="collision"):
        df.chem.to_mae(tmp_path / "ci.mae", index=True)


def test_secondary_molecule_full_key_passthrough(tmp_path):
    df = _mol_frame([_make_mol("primary")])
    # A full-key-named secondary column passes through verbatim (no double prefix).
    df["s_user_ligand"] = pd.Series(MoleculeArray([_make_mol("lig")]), dtype=MoleculeDtype())
    out = tmp_path / "secfk.mae"
    df.chem.to_mae(out)

    back = read_mae(out, tags=TAG_ALL)
    assert "s_user_ligand" in back.columns
    assert "s_user_s_user_ligand" not in back.columns


def test_secondary_collision_with_data_column_raises(tmp_path):
    df = _mol_frame([_make_mol("primary")])
    df["Other"] = pd.Series(MoleculeArray([_make_mol("s")]), dtype=MoleculeDtype())  # -> s_user_Other
    df["s_user_Other"] = pd.array(["x"], dtype="object")  # full key -> s_user_Other -> collision
    with pytest.raises(ValueError, match="collision"):
        df.chem.to_mae(tmp_path / "sc.mae")


def test_overwrite_through_symlink_updates_target_not_link(tmp_path):
    import os
    real = tmp_path / "real.mae"
    link = tmp_path / "link.mae"
    _mol_frame([_make_mol("orig")]).chem.to_mae(real)
    os.symlink(real, link)
    _mol_frame([_make_mol("updated")]).chem.to_mae(link)  # write through the symlink
    assert link.is_symlink(), "symlink must be preserved, not replaced by a regular file"
    assert os.path.realpath(link) == os.path.realpath(real)
    from oepandas_mae.reader import read_mae
    assert list(read_mae(real)["Title"]) == ["updated"]  # referent updated


def test_symlink_compression_mismatch_raises(tmp_path):
    import os
    real = tmp_path / "real.mae.gz"
    _mol_frame([_make_mol("x")]).chem.to_mae(real)   # real gz target
    link = tmp_path / "link.mae"                       # .mae link -> .mae.gz target
    os.symlink(real, link)
    with pytest.raises(ValueError, match="compression"):
        _mol_frame([_make_mol("y")]).chem.to_mae(link)


def test_symlink_matching_compression_ok(tmp_path):
    import os
    real = tmp_path / "real.mae.gz"
    _mol_frame([_make_mol("x")]).chem.to_mae(real)
    link = tmp_path / "link.mae.gz"                    # matching gz extension
    os.symlink(real, link)
    _mol_frame([_make_mol("y")]).chem.to_mae(link)     # must NOT raise
    assert link.is_symlink()
    from oepandas_mae.reader import read_mae
    assert list(read_mae(real)["Title"]) == ["y"]
