## dunders
__version__ = "1.0"
__author__ = "Isaac Overcast"

import os as _os
import atexit as _atexit

## Force matplotlib to behave on headless environments
import matplotlib
matplotlib.use("agg")

from . import gimmeSAD
