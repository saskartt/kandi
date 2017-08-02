#!/usr/bin/env python
import numpy as np

# Skip data from the beginning of the simulation before time-averaging
skip_time_avg = 2700.

# Profile heights
pr_heights_table = [1,2,4,6,8,10,12,14,16,17,18,20,25,30,100,128]


pr_heights_plot = np.arange(0,129,1)

# Building height
bld_height = 17.0
