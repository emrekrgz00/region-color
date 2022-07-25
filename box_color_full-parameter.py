
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

region_list = []

for box in patch_list:

    box_poly = Polygon(box.get_xy())
    region_list.append(box_poly)



cmap_plasma = py.cm.plasma

# fig = py.figure(figsize=(10, 10))

# ax = fig.add_subplot(1, 1, 1 , projection=wcs, facecolor="blue")

# grey_map = cm.get_cmap("gray").copy()
# grey_map.set_bad('white', 0.)
# im = ax.imshow(np.log10(hdu.data), cmap=grey_map, origin="lower", interpolation="bilinear")


# # for p in patch_list:
# #     ax.add_patch(p)
# # for t in artist_list:
# #     ax.add_artist(t)

# collection = PatchCollection(region_list, cmap=cmap_plasma, alpha=0.75)

# collection.set_array(colors)
# collection.set_clim(vmin=0.2, vmax=1.2)
# ax.add_collection(collection)


# ax.set_xlabel("Ra", fontsize=20)  #, weight=bold
# ax.set_ylabel("Dec", fontsize=20, labelpad=-2)  #, weight=bold
# ax.set_xlim(3900, 4300)
# ax.set_ylim(3900, 4300)
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.tick_params(which='both', width=1, labelsize=18, rasterized=True)
# ax.tick_params(which='major', length=1)
# ax.tick_params(which='minor', length=1)

# colorbar_ax = fig.colorbar(collection, orientation="horizontal", ax=ax, shrink=0.68)
# colorbar_ax.ax.set_xlabel(r"${\rm k T\,(keV)}$", size=19)
# colorbar_ax.ax.tick_params(which='both', width=1, labelsize=17)

# # divider = make_axes_locatable(ax)
# # colorbar_ax = divider.append_axes("top", "5%", pad="7%")

# # fig.colorbar(collection, cax=colorbar_ax,
# #             orientation="horizontal", format=r'%.2f')
# # colorbar_ax.set_xlabel(r"${\rm k T\,(keV)}$")
# # colorbar_ax.set_ylabel([])
# # colorbar_ax.set_yticks([])

# # fig.tight_layout()
# fig.savefig(cikti + "snr2_region_box_colored.png", dpi=300)
# py.close()








# logaritmik set
label_list = [r"${\rm k T\,(keV)}$", r"${\rm \log ((EM/10^{57})"+u"/Area)"+" \, (cm^{-3}\, arcsec^{-2})}$",
              r"${\rm \log (n_et) \,(cm^{-3}\,s)}$", r"${\rm O}$", r"${\rm Ne}$",
              r"${\rm Mg}$", r"${\rm Si}$", r"${\rm Fe}$"]

text_list = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]

vmin_list = [0.25, -0.8, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0]

vmax_list = [1.1, 0.4, 13.0, 0.4, 0.75, 0.5, 0.65, 1.7]

color_list = [kT, np.log10(em), np.log10(net), oxygen,
              neon, magnesium, silisium,
              fermium]


grey_map = cm.get_cmap("gray").copy()
grey_map.set_bad('white', 0.)

cmap = cm.get_cmap("jet").copy()
cmap.set_bad('blue', 1.)

cmap_plasma = cm.get_cmap("plasma").copy()

fig = py.figure(figsize=(17, 22))
for i in range(len(color_list)):

    ax_num = i + 1
    ax = fig.add_subplot(4, 2, ax_num, projection=wcs, facecolor="black")
    
    im = ax.imshow(np.log10(hdu.data), cmap=grey_map,
                   origin="lower", interpolation="bilinear")

    # for pa in patch_list:
    #     copy_pa = copy(pa)
    #     ax.add_patch(copy_pa)

    at = AnchoredText(text_list[i],
                      prop=dict(size=20),
                      frameon=True,
                      loc=1)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.1")
    ax.add_artist(at)

    ax.set_xlabel("Ra", fontsize=20)  #, weight=bold
    ax.set_ylabel("Dec", fontsize=20, labelpad=-2)  #, weight=bold
    ax.set_xlim(3900, 4300)
    ax.set_ylim(3900, 4300)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1, labelsize=18, rasterized=True)
    ax.tick_params(which='major', length=1)
    ax.tick_params(which='minor', length=1)

    colors = color_list[i]
    collection = PatchCollection(region_list, cmap=cmap_plasma, alpha=0.75)

    collection.set_array(colors)
    collection.set_clim(vmin=vmin_list[i], vmax=vmax_list[i])
    ax.add_collection(collection) 


    # # colorbar_ax = fig.add_axes([0.12, 0.965, 0.78, 0.02])

    colorbar_ax = fig.colorbar(collection, orientation="horizontal", ax=ax, shrink=0.45)
    # , ticks=np.arange(0, 2.75, 0.25))
    colorbar_ax.ax.set_xlabel(label_list[i], size=19) ## Alt label
    colorbar_ax.ax.tick_params(which='both', width=1, labelsize=17)


fig.canvas.draw()

# fig.subplots_adjust(left=-0.05, bottom=0.05, right=1.05, top=0.95, wspace=-0.20, hspace=0.1)

fig.tight_layout()
fig_name = cikti + "full_parameter-N132-region.png"
fig.savefig(fig_name, dpi=300)

py.close()