
from email import header
import os
import time
import numpy as np
import pylab as py
from astropy.coordinates import SkyCoord
from astropy.io import fits as pf
from astropy.io import ascii
from astropy import units as unit
from astropy.table import Table as ast_Table, Column, vstack, hstack
# from astroML.stats import binned_statistic_2d

import multiprocessing

from matplotlib.ticker import AutoMinorLocator

import glob

import pyregion
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle, Wedge
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import matplotlib.patheffects as pe

from regions import PixelRegion, PixCoord, PolygonPixelRegion


import numpy as np
from scipy.spatial import distance
from copy import copy

def sort_clockwise(points):

    # finding center points of polygon

    x_center = np.mean(points[:, 0])
    y_center = np.mean(points[:, 1])

    p_center = np.array([x_center, y_center]).reshape(1,2)

    dif_points = points - p_center

    degrees = np.rad2deg(np.arctan2(dif_points[:, 1], dif_points[:, 0]))

    sort_degrees = np.argsort(degrees)

    return points[sort_degrees]


t_start = time.time()

#ust klasor yapisi
data_yolu = os.pardir+"/data/"
cikti = os.pardir+"/cikti/"
datacikti = os.pardir+"/data_cikti/"
rapor = os.pardir+"/rapor/"
kod = os.getcwd()


hdu = pf.open(data_yolu + "snr2.fits")[0]
wcs = WCS(hdu.header, fix=False)

reg_file = data_yolu + "full-fk5.reg"

reg = pyregion.open(reg_file).as_imagecoord(hdu.header)

patch_list, artist_list = reg.get_mpl_patches_texts()

# kT,EM,net,O,Ne,Mg
data = ascii.read(data_yolu + "results_2-tek.csv", delimiter=",", guess=False)
kT = np.array(data["kT"])
em = np.array(data["EM"])
net = np.array(data["net"])
oxygen = np.array(data["O"])
neon = np.array(data["Ne"])
magnesium = np.array(data["Mg"])
silisium = np.array(data["Si"])
fermium = np.array(data["Fe"])

colors = kT



cmap_plasma = py.cm.plasma

fig = py.figure(figsize=(10, 10))

ax = fig.add_subplot(1, 1, 1 , projection=wcs, facecolor="blue")

grey_map = cm.gray
grey_map.set_bad('white', 0.)
im = ax.imshow(np.log10(hdu.data), cmap=grey_map,
               origin="lower", interpolation="bilinear")


for p in patch_list:
    ax.add_patch(p)
for t in artist_list:
    ax.add_artist(t)

ax.set_xlabel("Ra", fontsize=20)  #, weight=bold
ax.set_ylabel("Dec", fontsize=20)  #, weight=bold
ax.set_xlim(3900, 4300)
ax.set_ylim(3900, 4300)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=1, labelsize=18, rasterized=True)
ax.tick_params(which='major', length=1)
ax.tick_params(which='minor', length=1)

# fig.tight_layout()
fig.savefig(cikti + "snr2_region_box.png", dpi=300)
py.close()