import pandas as pd
from oepandas.arrays import MoleculeDtype

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
