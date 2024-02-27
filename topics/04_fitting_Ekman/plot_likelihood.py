import numpy as np
import xarray as xr

import traceback
from pathlib import Path

import argparse



parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-file', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--x-rng', type=float, nargs=2, help='The x axis range to be plot in km.', default=[None, None])
parser.add_argument('--no-display', action="store_true")

args = parser.parse_args()
print(args)

ds = xr.open_dataset(args.input_file)

coords = { varname : ds.coords[varname].to_numpy() for varname in ds.coords}

# Plot data
print("Loading Matplotlib...")
import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    mpl.rc('font', size=15)
    mpl.rc('axes', labelsize=15)

print("Done.")     
 
 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import tool_fig_config

ncol = 1
nrow = 2

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 6,
    h = [3, 4],
    wspace = 1.0,
    hspace = .7,
    w_left = 1.0,
    w_right = 1.5,
    h_bottom = 1.0,
    h_top = 1.0,
    ncol = ncol,
    nrow = nrow,
)


fig, ax = plt.subplots(
    nrow, ncol,
    figsize=figsize,
    subplot_kw=dict(aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
)


_ax = ax.flatten()[0]
_ax.plot(ds.coords["E0"], ds["loglikelihood"], "k-", linewidth=2)

trans = transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
_ax.plot([ds.attrs["E0_maxlikelihood"]]*2, [0, 1], "r-", linewidth=1, transform=trans)


_ax.set_title("log likelihood")
_ax.set_xlabel("$ E_0 $ [ $\\mathrm{s}^{-1}$ ]")

_ax = ax.flatten()[1]

for i in range(len(ds.coords["s_T"])):
    _ax.plot([ds["u_0"][i], ds["obs_U"][i]], [ds["v_0"][i], ds["obs_V"][i]], color="gray", alpha=0.8, linewidth=1.0)


_ax.plot(ds["u_0"], ds["v_0"],     "r--", linewidth=2, label="fitted", marker = '.', markersize = 10)
_ax.plot(ds["obs_U"], ds["obs_V"], "k-", linewidth=2, label="WRF", marker = '.', markersize = 10)

_ax.legend()

_ax.set_xlabel("$ U $ [ $ \\mathrm{m} / \\mathrm{s}$ ]")
_ax.set_ylabel("$ V $ [ $ \\mathrm{m} / \\mathrm{s}$ ]")

for _ax in ax.flatten():
    _ax.grid()

if args.output != "":
    print("Saving output: ", args.output)    
    fig.savefig(args.output, dpi=200)

if not args.no_display:
    plt.show()

