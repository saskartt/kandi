#!/usr/bin/env python


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

'''
Generates a simple discretized colorbar and saves it in vector format.
'''

fig = plt.figure(figsize=(1, 4))
ax = fig.add_axes([0.02, 0.025, 0.4, 0.85])

# The third example illustrates the use of custom length colorbar
# extensions, used on a colorbar with discrete intervals.
cmap = mpl.cm.get_cmap("RdBu_r")
cmap.set_over((1., 0., 0.))
cmap.set_under((0., 0., 1.))
cmap.set_bad('black',1.)

bounds = np.linspace(-0.75,0.75,16)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                boundaries=bounds,
                                # Make the length of each extension
                                # the same as the length of the
                                # interior colors:2
                                extendfrac='auto',
                                ticks=bounds,
                                spacing='uniform',
                                orientation='vertical')
ax.set_title('$\Delta |\mathbf{w}| \; \mathrm{[m/s]}$')
cb.ax.tick_params(labelsize=10)

fig.savefig("colorbar.svg")
