#!/usr/bin/env python
import sys
import argparse
import numpy as np
from PIL import Image
from kandiLib import *
from mapTools import *

'''
Generates a Numpy Z raster file from TIFF image data
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='createCanopyRaster.py', description='''Calcluates various quantities from profiles.''')
parser.add_argument("-f", "--file", type=str, help="Name of the input TIFF file.")
parser.add_argument("-fo", "--fileout", type=str, help="Name of the output raster data file.")
parser.add_argument("-ch", "--canopyheight", type=float, help="Canopy height.")
args = parser.parse_args()

#==========================================================#

# try
img = Image.open(args.file)
# except:
  # raise IOError("Input file {} not found!".format(args.file))

imgarr = np.array(img)

imgarr[np.nonzero(imgarr)] = args.canopyheight

Rdict={'R' : imgarr, 'GlobOrig' : [0.0,0.0], 'dPx' : [1.0, 1.0]}
saveTileAsNumpyZ(args.fileout, Rdict)
