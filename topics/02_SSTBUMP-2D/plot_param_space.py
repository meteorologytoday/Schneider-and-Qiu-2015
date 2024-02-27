import numpy as np
import xarray as xr

import traceback
from pathlib import Path

import argparse



parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input directory', required=True)
parser.add_argument('--selected-Ug', type=float, help='Selected Ug', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--no-display', action="store_true")

args = parser.parse_args()
print(args)

pathlist = []
sorting = []
for i, path in enumerate(Path(args.input_dir).glob('*.nc')):

    _ds = xr.open_dataset(path)
    if _ds.attrs["Ug"] != args.selected_Ug:
        continue


    pathlist.append((path, _ds.attrs["dSST"],))
    #sorting.append(_ds.attrs["dSST"])

pathlist.sort(key=lambda x: x[1])
pathlist = [ path[0] for path in pathlist]

print("Going to load following files: ")
for i, path in enumerate(pathlist):
    print("[%d] %s" % ( i, path, ))

N = len(pathlist)
varnames = ["dSST", "U_mean", "delta_mean", "Udelta_CORR", "Udelta_FULL"]

data = { varname : np.zeros((N,)) for varname in varnames }

for i, path in enumerate(pathlist):
    
    print("Processing file: %s. Ug = %.1f, dSST = %.1f" % (path, _ds.attrs["Ug"], _ds.attrs["dSST"],))

    _ds = xr.open_dataset(path)
    _ds = _ds.isel(y_T=0)
    
    data["dSST"][i] = _ds.attrs["dSST"]

    # Do integration
    dx_T = _ds["dx_T"].to_numpy()
   
    U = ( (_ds["u_0"] + _ds["u_1"])**2 + (_ds["v_0"] + _ds["v_1"])**2 )**0.5
    U = U.isel(z_T=0).to_numpy() 
    delta = (_ds["sst_1"] -  _ds["pt_1"]).to_numpy()
 
    #U_1 = ((_ds["u_1"]**2 + _ds["v_1"]**2)**0.5).isel(z_T=0).to_numpy()
    #U_1 = _ds["u_1"].isel(z_T=0).to_numpy()
    #U_0 = U_1*0 + (( _ds["u_0"]**2 + _ds["v_0"]**2 )**0.5).isel(z_T=0).to_numpy()

    #data["U1delta1"][i] = np.sum(dx_T * U_1 * delta_1) / np.sum(dx_T)
    #data["U0delta1"][i] = np.sum(dx_T * U_0 * delta_1) / np.sum(dx_T)
    U_mean = np.sum(dx_T * U) / np.sum(dx_T)
    delta_mean = np.sum(dx_T * delta) / np.sum(dx_T)

    data["Udelta_FULL"][i] = np.sum(dx_T * U * delta) / np.sum(dx_T)
    data["U_mean"][i] = U_mean
    data["delta_mean"][i] = delta_mean
    data["Udelta_CORR"][i] = np.sum(dx_T * (U - U_mean) * (delta - delta_mean)) / np.sum(dx_T)
    

const_H = 6.9
ref_idx = 0
data_dif = dict(
    WND = np.zeros((N,)),
    THM = np.zeros((N,)),
    COR = np.zeros((N,)),
    FUL = np.zeros((N,)),
)

for i, path in enumerate(pathlist):
    
    data_dif["WND"][i] = const_H * (data["U_mean"][i] - data["U_mean"][ref_idx]) * data["delta_mean"][ref_idx]
    data_dif["THM"][i] = const_H * (data["delta_mean"][i] - data["delta_mean"][ref_idx]) * data["U_mean"][ref_idx]
    data_dif["COR"][i] = const_H * (data["Udelta_CORR"][i] - data["Udelta_CORR"][ref_idx])
    data_dif["FUL"][i] = const_H * (data["Udelta_FULL"][i] - data["Udelta_FULL"][ref_idx])

data_dif["RES"] = data_dif["FUL"] - (data_dif["THM"] + data_dif["COR"] + data_dif["WND"])

print("Domain size: sum(dx_T) = ", np.sum(dx_T))


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

varnames = ["WND", "THM", "COR", "RES"]


ncol = 1
nrow = len(varnames)

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 6,
    h = [3,] * nrow,
    wspace = 1.0,
    hspace = 1.0,
    w_left = 1.0,
    w_right = 2.0,
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
    sharex=True,
)

for i, varname in enumerate(varnames):

    _ax = ax[i, 0]
    _ax.plot(data["dSST"], data_dif[varname])
    _ax.set_title("Variable: %s" % (varname,))

    _ax.set_xlabel("$ \\Delta \\mathrm{SST}$ [ K ]")

for _ax in ax.flatten():
    _ax.grid()

if not args.no_display:
    plt.show()

if args.output != "":

    print("Saving output: ", args.output)    
    fig.savefig(args.output, dpi=200)


