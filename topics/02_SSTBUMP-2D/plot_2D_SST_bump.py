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


ds = ds.isel(y_T=0)
#print(ds)




coords = { varname : ds.coords[varname].to_numpy() for varname in ds.coords}
coords["x_T"] /= 1e3

#coords["pz_T"] = coords["z_T"] * ds.attrs["h_0"]
coords["pz_T"] = coords["z_T"][:, None] * ( ds.h_1.to_numpy()[None, :] + ds.attrs["h_0"] )
coords["pz_W"] = coords["z_W"][:, None] * ( ds.h_1.to_numpy()[None, :] + ds.attrs["h_0"] )

x_TT = np.repeat(np.reshape(coords["x_T"], (1, -1)), len(coords["z_T"]), axis=0)
x_TW = np.repeat(np.reshape(coords["x_T"], (1, -1)), len(coords["z_W"]), axis=0)

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
nrow = 3

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 6,
    h = [3, 2, 2],
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
    sharex=True,
)


w_1_levels = np.linspace(-1, 1, 11) * 0.5

w_1 = ds.w_1 * 1e2 * ds.attrs["h_0"]
mappable = ax[0,0].contourf(x_TW, coords["pz_W"], w_1, w_1_levels, cmap="bwr", extend="both")
cs = ax[0,0].contour(x_TW, coords["pz_W"], w_1, w_1_levels, colors="black", linestyles="solid", linewidths=0.5)
cax = tool_fig_config.addAxesNextToAxes(fig, ax[0, 0], "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.00)

ax[0,0].plot(coords["x_T"], ds.h_1 + ds.attrs["h_0"], "k--", label="$h^{\\left(1\\right)}$")

#plot_info = plot_infos[varnames[i]]
#unit_str = "" if plot_info["unit"] == "" else " [ %s ]" % (plot_info["unit"],)
cb.ax.set_ylabel("$ w^{\\left(1\\right)}$ [ cm / s ]", size=15)

ax[1,0].plot(coords["x_T"], ds.u_1.isel(z_T=0), color="black", label="$ u^{\\left(1\\right)}_\\mathrm{sfc}$")
ax[1,0].plot(coords["x_T"], ds.v_1.isel(z_T=0), color="red", label="$ v^{\\left(1\\right)}_\\mathrm{sfc}$")

ax[1,0].set_title("$\\left( u^{\\left(0\\right)}_\\mathrm{sfc}, v^{\\left(0\\right)}_\\mathrm{sfc} \\right) = \\left( %.2f, %.2f \\right)$" % (ds.u_0.isel(z_T=0), ds.v_0.isel(z_T=0)))

sfc_u_1 = ds.u_1.isel(z_T=0).to_numpy()
sfc_v_1 = ds.v_1.isel(z_T=0).to_numpy()
sfc_u_0 = ds.u_0.isel(z_T=0).to_numpy()
sfc_v_0 = ds.v_0.isel(z_T=0).to_numpy()
VEL_tot_w_SST  = np.sqrt((sfc_u_0 + sfc_u_1)**2 + (sfc_v_0 + sfc_v_1)**2)
VEL_tot_wo_SST  = np.sqrt((sfc_u_0 + sfc_u_1*0)**2 + (sfc_v_0 + sfc_v_1*0)**2)

#ax[2,0].plot(coords["x_T"], VEL_tot_w_SST,  "k-",  label="total $u$ w/ SST bump")
#ax[2,0].plot(coords["x_T"], VEL_tot_wo_SST, "k--", label="total $u$ wo/ SST bump")

#ax[2,0].plot(coords["x_T"], ds.h_1, "k-", label="$h^{\\left(1\\right)}$")


ax[2,0].plot(coords["x_T"], ds.sst_1, color="blue", label="$ \\mathrm{SST}^{\\left(1\\right)}$")
ax[2,0].plot(coords["x_T"], ds.pt_1, color="red", label="$ \\Theta^{\\left(1\\right)}$")

ax[1,0].legend(ncol=2)
ax[2,0].legend()
#ax[3,0].legend()

ax[0,0].set_ylabel("z [ m ]")
ax[1,0].set_ylabel("[ m / s ]")
#ax[2,0].set_ylabel("[ m ]")
ax[2,0].set_ylabel("[ K ]")

ax[0,0].set_ylim([0, ds.attrs["h_0"] * 1.1])

ax[0, 0].set_xlim(args.x_rng)
for _ax in ax.flatten():
    _ax.grid()
    _ax.set_xlabel("x [ km ]")

if args.output != "":

    print("Saving output: ", args.output)    
    fig.savefig(args.output, dpi=200)

if not args.no_display:
    plt.show()

