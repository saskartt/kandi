#!/usr/bin/env python
import numpy as np

# Skip data from the beginning of the simulation before time-averaging
timespan = [1800, 7200.]

# Profile heights
pr_heights_table = [1,2,4,6,8,10]


pr_heights_plot = np.arange(0,129,1)

# Building height
bld_height = 16.0

# Canopy settings
tree_location = [27.0, 37.0, 91.0, 101.0]
tree_height = [11.0, 11.0, 11.0, 11.0]
tree_radius = [4.0, 4.0, 4.0, 4.0]

# Reference friction velocity
z_scale = 100 # Characteristic length scale in meters
p_grad = -0.003# Pressure gradient
rho = 1.25 # Air density
u_ref=np.sqrt(-(1./rho)*p_grad*z_scale)
