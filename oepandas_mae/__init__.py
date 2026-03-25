import logging

__version__ = "0.1.0"

# ---- Logging ----

log = logging.getLogger("oepandas-mae")


class LevelSpecificFormatter(logging.Formatter):
    """A logging formatter that prefixes every level but INFO."""

    NORMAL_FORMAT = "%(message)s"
    LEVEL_SPECIFIC_FORMAT = "%(levelname)s: %(message)s"

    def __init__(self):
        super().__init__(fmt=self.NORMAL_FORMAT, datefmt=None, style="%")

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno == logging.INFO:
            self._style._fmt = self.NORMAL_FORMAT
        else:
            self._style._fmt = self.LEVEL_SPECIFIC_FORMAT
        return logging.Formatter.format(self, record)


ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(LevelSpecificFormatter())
log.addHandler(ch)
log.setLevel(logging.INFO)

# ---- Public API ----

from .reader import read_mae  # noqa: E402
from ._compat import (  # noqa: E402
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
)

# ---- Auto-registration ----

import pandas as pd  # noqa: E402
pd.read_mae = read_mae

try:
    import oepandas  # noqa: E402
    oepandas.read_mae = read_mae
except ImportError:
    pass

# ---- Exports ----

__all__ = [
    "__version__",
    "read_mae",
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
]
