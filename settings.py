#!/usr/bin/env python
import numpy as np

# Skip data from the beginning of the simulation before time-averaging
timespan = [3600., 7200.]

# Profile heights
pr_heights_table = [1,2,4,6,8,10,12,14,16,17,18,20,25,30,100,128]


pr_heights_plot = np.arange(0,129,1)

# Building height
bld_height = 16.0

# Canopy settings
tree_location = [27.0, 37.0]
tree_height = [11.0, 11.0]
tree_radius = [4.0, 4.0]
