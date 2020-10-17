###############################################################################
# linkFiles.py
# 
# Last update: 6/18/2020 by Isaac Arseneau
#
# Description: Links the necessary files to recreate the wrf runs
# along with the traj files generated on that day
###############################################################################

### Imports ###

from netCDF4 import Dataset
import numpy as np
import trajUtils_v as utils
from metpy.units import units
import sys
from os import (symlink, unlink, path)
args = sys.argv


### Setup/Parameters ###

utils.getRecentDate()