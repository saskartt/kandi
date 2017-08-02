#!/usr/bin/env python
import sys
import argparse
import numpy as np
from PIL import Image
import matplotlib as mpl
# mpl.use('GTKCairo')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from settings import *

'''
Plot raster image with custom colorbar and axis.
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='createCanopyRaster.py', description='''Calcluates various quantities from profiles.''')
parser.add_argument("-f", "--file", type=str, help="Name of the input TIFF file.")
parser.add_argument("-s", "--save", type=str, help="Save resulting figure as.")
parser.add_argument("-cm", "--colormap", type=str, help="Colormap style")
parser.add_argument("-l", "--lims", type=float, nargs=2, help="Colorbar limits.")
parser.add_argument("-y", "--ymax", type=float, help="Maximum (relative) height.")
parser.add_argument("-hl", "--height", type=float, help="Building height in pixels.")
args = parser.parse_args()

#==========================================================#

try:
  img = Image.open(args.file)
except:
  raise IOError("Input file {} not found!".format(args.file))

imgarr = np.array(img)
if (args.colormap):
  cmap = mpl.cm.get_cmap(args.colormap)
  norm = mpl.colors.Normalize(vmin=args.lims[0], vmax=args.lims[1])
else:
  cmap = None
  norm = None

hlen = args.height
xlen = np.shape(imgarr)[1]
ylen = np.shape(imgarr)[0]
xticklabels = np.arange(-1.5,2.0,0.5)
xticks = np.arange(xlen/2-1.5*hlen,xlen/2+2*hlen,hlen/2)

print(ylen)

yticks = np.arange(0.0,3.0*hlen+hlen/2,hlen/2)
yticklabels = np.arange(0.0,3.0,0.5)
print(np.shape(imgarr))

fig=plt.figure()
ax = plt.gca()

imgplot = plt.imshow(imgarr,cmap=cmap, norm=norm,extent=[0,xlen,0,ylen],aspect=1)
plt.xticks(xticks, xticklabels)
plt.yticks(yticks, yticklabels)
ax.set_ylim([0,args.ymax*hlen])
# ax.invert_yaxis()
 ###################

if (cmap):
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cbar = fig.colorbar(imgplot,cax=cax)

####################
if (args.save):
  fig.savefig( args.save, format='eps', dpi=600, bbox_inches='tight')
  print("Figure {} saved.".format(args.save))
plt.show()
