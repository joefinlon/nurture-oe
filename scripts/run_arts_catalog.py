"""Task (c): ARTS SSDB v1.1 catalog (``SSD/TotallyRandom``), incl. Ka-band frequencies.

Writes ``docs/reports/arts_catalog.csv`` (one row per habit) and
``docs/reports/arts_catalog_sizes.csv`` (one row per size file). Read-only on the data store.
"""

import sys

from nurture_oe.config import REPO_ROOT, load_paths
from nurture_oe.ssp.arts_catalog import build_catalog


def main():
    cat, sizes = build_catalog(load_paths()["arts_ssdb"])
    rep = REPO_ROOT / "docs" / "reports"
    cat.to_csv(rep / "arts_catalog.csv", index=False, float_format="%.6g")
    sizes.to_csv(rep / "arts_catalog_sizes.csv", index=False, float_format="%.6g")
    print(cat.to_string())


if __name__ == "__main__":
    sys.exit(main())
