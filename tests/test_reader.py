import pandas as pd
import pytest
from openeye import oechem
from oepandas.arrays import MoleculeDtype

from oepandas_mae._compat import (
    TAG_ALL,
    TAG_NAME,
    TAG_NONE,
    PERCEPTION_NONE,
    OEMaestroReaderConfig,
)
from oepandas_mae.reader import read_mae


class TestReadMaeBasic:
    def test_basic_read(self, simple_mae):
        df = read_mae(simple_mae)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "Molecule" in df.columns
        assert "Title" in df.columns
        assert isinstance(df["Molecule"].dtype, MoleculeDtype)
        assert pd.api.types.is_string_dtype(df["Title"].dtype)

    def test_basic_read_has_data_columns(self, protein_maegz):
        df = read_mae(protein_maegz)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        # Should have columns beyond just Molecule and Title (8G66.maegz has CT properties)
        assert len(df.columns) > 2

    def test_empty_file_returns_empty_dataframe(self, tmp_path):
        empty_mae = tmp_path / "empty.mae"
        empty_mae.write_text("")
        df = read_mae(empty_mae)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert "Molecule" in df.columns
        assert "Title" in df.columns


class TestReadMaeConfig:
    def test_default_tags_is_tag_name(self, protein_maegz):
        """Default tag format should be TAG_NAME (clean column names without type/owner prefix)."""
        df = read_mae(protein_maegz)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0, "Need data columns to test tag formatting"
        for col in data_cols:
            # TAG_ALL columns look like "r_m_something" or "s_m_something"
            assert not (len(col) > 4 and col[1] == "_" and col[3] == "_"), \
                f"Column '{col}' looks like TAG_ALL format, expected TAG_NAME"

    def test_tags_all_produces_full_keys(self, protein_maegz):
        """TAG_ALL should produce full Maestro-style keys like r_m_pdb_tfactor."""
        df = read_mae(protein_maegz, tags=TAG_ALL)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0
        # At least some columns should have the t_o_ prefix pattern
        has_full_key = any(
            len(c) > 4 and c[1] == "_" and c[3] == "_" for c in data_cols
        )
        assert has_full_key, f"Expected TAG_ALL columns, got: {data_cols}"

    def test_tags_none_produces_no_data_columns(self, protein_maegz):
        """TAG_NONE should produce no data columns."""
        df = read_mae(protein_maegz, tags=TAG_NONE)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) == 0

    def test_perception_none(self, simple_mae):
        """PERCEPTION_NONE should still read molecules successfully."""
        df = read_mae(simple_mae, perception=PERCEPTION_NONE)
        assert len(df) > 0
        assert "Molecule" in df.columns

    def test_config_object(self, protein_maegz):
        """Config object should be accepted."""
        cfg = OEMaestroReaderConfig()
        cfg.tags = TAG_ALL
        df = read_mae(protein_maegz, config=cfg)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        has_full_key = any(
            len(c) > 4 and c[1] == "_" and c[3] == "_" for c in data_cols
        )
        assert has_full_key

    def test_kwargs_override_config(self, protein_maegz):
        """Kwargs should override config values."""
        cfg = OEMaestroReaderConfig()
        cfg.tags = TAG_ALL
        # tags kwarg should override config
        df = read_mae(protein_maegz, config=cfg, tags=TAG_NAME)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0, "Need data columns to test override"
        for col in data_cols:
            assert not (len(col) > 4 and col[1] == "_" and col[3] == "_"), \
                f"Column '{col}' looks like TAG_ALL, but TAG_NAME kwarg should have overridden config"


class TestReadMaeColumns:
    def test_custom_molecule_column_name(self, simple_mae):
        df = read_mae(simple_mae, molecule_column="Mol")
        assert "Mol" in df.columns
        assert "Molecule" not in df.columns
        assert isinstance(df["Mol"].dtype, MoleculeDtype)

    def test_custom_title_column_name(self, simple_mae):
        df = read_mae(simple_mae, title_column="Name")
        assert "Name" in df.columns
        assert "Title" not in df.columns

    def test_no_title_suppresses_title_column(self, simple_mae):
        df = read_mae(simple_mae, no_title=True)
        assert "Title" not in df.columns
        assert "Molecule" in df.columns

    def test_title_column_none_suppresses_title(self, simple_mae):
        df = read_mae(simple_mae, title_column=None)
        assert "Title" not in df.columns
        assert "Molecule" in df.columns


class TestReadMaeGzip:
    def test_mae_gz(self, simple_mae, simple_mae_gz):
        """Reading .mae.gz should produce the same result as .mae."""
        df_plain = read_mae(simple_mae)
        df_gz = read_mae(simple_mae_gz)
        assert len(df_plain) == len(df_gz)
        assert list(df_plain.columns) == list(df_gz.columns)

    def test_maegz(self, protein_maegz):
        """Reading .maegz should work."""
        df = read_mae(protein_maegz)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "Molecule" in df.columns


class TestReadMaeUsecols:
    def test_usecols_filters_columns(self, protein_maegz):
        """usecols should limit which data columns appear."""
        df_all = read_mae(protein_maegz)
        data_cols = [c for c in df_all.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0, "Need data columns to test usecols"

        target_col = data_cols[0]
        df_filtered = read_mae(protein_maegz, usecols=[target_col])
        filtered_data_cols = [c for c in df_filtered.columns if c not in ("Molecule", "Title")]
        assert filtered_data_cols == [target_col]

    def test_usecols_string(self, protein_maegz):
        """usecols as a single string should work."""
        df_all = read_mae(protein_maegz)
        data_cols = [c for c in df_all.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0
        target_col = data_cols[0]
        df_filtered = read_mae(protein_maegz, usecols=target_col)
        filtered_data_cols = [c for c in df_filtered.columns if c not in ("Molecule", "Title")]
        assert target_col in filtered_data_cols


class TestReadMaeNumeric:
    def test_numeric_conversion(self, protein_maegz):
        """numeric should convert string columns to numeric types."""
        df = read_mae(protein_maegz, tags=TAG_NAME)
        data_cols = [c for c in df.columns if c not in ("Molecule", "Title")]
        assert len(data_cols) > 0, "Need data columns to test numeric conversion"
        # Find a numeric-looking column by checking if values can be converted
        numeric_col = None
        for col in data_cols:
            try:
                pd.to_numeric(df[col], errors="raise")
                numeric_col = col
                break
            except (ValueError, TypeError):
                continue

        assert numeric_col is not None, "Need at least one numeric-convertible column"
        df2 = read_mae(protein_maegz, tags=TAG_NAME, numeric=numeric_col)
        assert pd.api.types.is_numeric_dtype(df2[numeric_col])

    def test_numeric_molecule_column_raises(self, simple_mae):
        """Passing the molecule column to numeric should raise KeyError."""
        with pytest.raises(KeyError, match="Cannot make molecule column"):
            read_mae(simple_mae, numeric="Molecule")


class TestReadMaeSmiles:
    def test_add_smiles_true(self, simple_mae):
        """add_smiles=True should add a SMILES column."""
        df = read_mae(simple_mae, add_smiles=True)
        assert "Molecule SMILES" in df.columns
        assert df["Molecule SMILES"].dtype == object
        assert all(isinstance(s, str) for s in df["Molecule SMILES"] if pd.notna(s))

    def test_add_smiles_false(self, simple_mae):
        """add_smiles=False should not add a SMILES column."""
        df = read_mae(simple_mae, add_smiles=False)
        smiles_cols = [c for c in df.columns if "SMILES" in c]
        assert len(smiles_cols) == 0

    def test_add_smiles_none(self, simple_mae):
        """add_smiles=None (default) should not add a SMILES column."""
        df = read_mae(simple_mae)
        smiles_cols = [c for c in df.columns if "SMILES" in c]
        assert len(smiles_cols) == 0


class TestReadMaeConformerGrouping:
    def test_absolute_groups_conformers(self, conformers_mae):
        """Absolute conformer test groups consecutive identical molecules."""
        df = read_mae(conformers_mae, conformer_test="absolute")
        # conformers.mae has 2 CTs both titled "Ethanol" with identical connectivity
        assert len(df) == 1
        mol = df["Molecule"].iloc[0]
        assert mol.NumConfs() == 2

    def test_different_titles_not_grouped(self, multi_mae):
        """Molecules with different titles should not be grouped even with default conf test."""
        df = read_mae(multi_mae)
        for mol in df["Molecule"]:
            if mol is not None:
                assert isinstance(mol, oechem.OEMol)
                assert mol.NumConfs() == 1

    def test_invalid_conformer_test_raises(self, simple_mae):
        with pytest.raises(ValueError, match="Invalid conformer_test"):
            read_mae(simple_mae, conformer_test="invalid")

    def test_all_conformer_test_values_accepted(self, simple_mae):
        """All valid conformer_test values should be accepted without error."""
        for ct in ("default", "absolute", "absolute_canonical", "isomeric", "omega"):
            df = read_mae(simple_mae, conformer_test=ct)
            assert isinstance(df, pd.DataFrame)


class TestReadMaeRegistration:
    def test_pd_read_mae_registered(self):
        """pd.read_mae should be available after importing oepandas_mae."""
        import oepandas_mae  # noqa: F401
        assert hasattr(pd, "read_mae")
        assert callable(pd.read_mae)

    def test_oepandas_read_mae_registered(self):
        """oepandas.read_mae should be available after importing oepandas_mae."""
        import oepandas
        import oepandas_mae  # noqa: F401
        assert hasattr(oepandas, "read_mae")
        assert callable(oepandas.read_mae)

    def test_registered_function_works(self, simple_mae):
        """The registered function should work identically to the direct import."""
        import oepandas_mae  # noqa: F401
        df = pd.read_mae(simple_mae)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
