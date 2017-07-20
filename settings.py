#!/usr/bin/env python
import numpy as np

# Skip data from the beginning of the simulation before time-averaging
skip_time_avg = 3600.

# Profile heights
pr_heights_table = [1,2,4,6,8,10,12,14,16,18,20,25,30,100,128]
pr_table_header = ['Variable','1 m', '2 m', '4 m', '6 m', '8 m' , '10 m', '12 m', '14 m', '16 m' , '18 m', '20 m', '25 m', '30 m', '100 m', '128 m']

pr_heights_plot = np.arange(0,129,1)
