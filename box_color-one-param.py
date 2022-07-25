
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
import matplotlib as pml
from matplotlib.ticker import AutoMinorLocator

import glob

import pyregion
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle, Wedge
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

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
data = ascii.read(data_yolu + "results_chi.csv", delimiter=",", guess=False)
kT = np.array(data["nH1"])
# em = np.array(data["EM"])
# net = np.array(data["net"])
# oxygen = np.array(data["O"])
# neon = np.array(data["Ne"])
# magnesium = np.array(data["Mg"])
# silisium = np.array(data["Si"])
# fermium = np.array(data["Fe"])

colors = kT

region_list = []

for box in patch_list:

    box_poly = Polygon(box.get_xy())
    region_list.append(box_poly)



cmap_plasma = py.cm.plasma

fig = py.figure(figsize=(10, 10))

ax = fig.add_subplot(1, 1, 1 , projection=wcs, facecolor="blue")

grey_map = cm.get_cmap("gray").copy()
grey_map.set_bad('white', 0.)
im = ax.imshow(np.log10(hdu.data), cmap=grey_map, origin="lower", interpolation="bilinear")
                # np.log10(hdu.data)

# for p in patch_list:
#     ax.add_patch(p)
# for t in artist_list:
#     ax.add_artist(t)

collection = PatchCollection(region_list, cmap=cmap_plasma, alpha=0.75)

collection.set_array(colors)
collection.set_clim(vmin=1.0, vmax=55.8)
ax.add_collection(collection)


ax.set_xlabel("Ra", fontsize=20)  #, weight=bold
ax.set_ylabel("Dec", fontsize=20, labelpad=-2)  #, weight=bold
ax.set_xlim(3900, 4300)
ax.set_ylim(3900, 4300)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=1, labelsize=18, rasterized=True)
ax.tick_params(which='major', length=1)
ax.tick_params(which='minor', length=1)

colorbar_ax = fig.colorbar(collection, orientation="horizontal", ax=ax, shrink=0.68)
colorbar_ax.ax.set_xlabel(r"${\rm k T\,(keV)}$", size=19)
colorbar_ax.ax.tick_params(which='both', width=1, labelsize=17)

# divider = make_axes_locatable(ax)
# colorbar_ax = divider.append_axes("top", "5%", pad="7%")

# fig.colorbar(collection, cax=colorbar_ax,
#             orientation="horizontal", format=r'%.2f')
# colorbar_ax.set_xlabel(r"${\rm k T\,(keV)}$")
# colorbar_ax.set_ylabel([])
# colorbar_ax.set_yticks([])

# fig.tight_layout()
fig.savefig(cikti + "snr2_region_box_colored-nh3.png", dpi=300)
py.close()