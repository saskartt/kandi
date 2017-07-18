#!/usr/bin/env python
import sys
import argparse
import numpy as np
import netCDF4 as nc
import tabulate as tb
from settings import *
from kandiLib import *

'''
Calculate various values from profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='analyzeProfiles.py', description='''Calcluates various quantities from profiles.''')
parser.add_argument("rfile", type=str, nargs='?', default=None, help="Name of the input NetCDF4 profile file.")
args = parser.parse_args()

#==========================================================#

# Read in the dataset
ds = nc.Dataset(args.rfile)

# Skip first hour of the simulation
timepoints, = np.where(ds.variables['time'][:]>=skip_time_avg)
heights = [1,2,4,6,8,10,12,14,16,18,20,25,30]
header = ['Variable','1 m', '2 m', '4 m', '6 m', '8 m' , '10 m', '12 m', '14 m', '16 m' , '18 m', '20 m', '25 m', '30 m']


# Calculate average profiles for wind component
sep = [[''] * (len(pr_heights)+1)] # empty line string
print("#=================== Statistical domain 00 ===================#")
pr_01 = averageProfilesWS('00', timepoints, pr_heights, ds)
pr_02 = averageProfilesVariances('00', timepoints, pr_heights,ds)
qty_list = np.concatenate((pr_01, sep, pr_02), axis=0)
print(tb.tabulate(qty_list, headers=header)+"\n")
print("#=================== Statistical domain 01 ===================#")
pr_01 = averageProfilesWS('01', timepoints, pr_heights, ds)
pr_02 = averageProfilesVariances('01', timepoints, pr_heights,ds)
qty_list = np.concatenate((pr_01, sep, pr_02), axis=0)
print(tb.tabulate(qty_list, headers=header)+"\n")
print("#=================== Statistical domain 02 ===================#")
pr_01 = averageProfilesWS('02', timepoints, pr_heights, ds)
pr_02 = averageProfilesVariances('02', timepoints, pr_heights,ds)
qty_list = np.concatenate((pr_01, sep, pr_02), axis=0)
print(tb.tabulate(qty_list, headers=header)+"\n")
