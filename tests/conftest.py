from pathlib import Path

import pytest

ASSETS_DIR = Path(__file__).parent / "assets"


@pytest.fixture
def simple_mae() -> Path:
    return ASSETS_DIR / "simple.mae"


@pytest.fixture
def simple_mae_gz() -> Path:
    return ASSETS_DIR / "simple.mae.gz"


@pytest.fixture
def multi_mae() -> Path:
    return ASSETS_DIR / "multi.mae"


@pytest.fixture
def protein_mae() -> Path:
    return ASSETS_DIR / "protein.mae"


@pytest.fixture
def conformers_mae() -> Path:
    return ASSETS_DIR / "conformers.mae"


@pytest.fixture
def protein_maegz() -> Path:
    return ASSETS_DIR / "8G66.maegz"
