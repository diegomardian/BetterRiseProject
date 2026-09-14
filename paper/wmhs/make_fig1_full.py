#!/usr/bin/env python
"""Render the full paper's candidate label from the existing calibration table.

Run from the repo root. No simulation is performed; the original figure used
by the extended abstract is preserved.
"""
from pathlib import Path

from make_fig1 import main

if __name__ == "__main__":
    raise SystemExit(main(out=Path("paper/wmhs/figures/full/fig1_calibration.pdf"),
                         annotation="candidate"))
