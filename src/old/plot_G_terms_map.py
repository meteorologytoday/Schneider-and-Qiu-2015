import numpy as np
import xarray as xr

import traceback
from pathlib import Path
import argparse
from datetime import (timedelta, datetime, timezone)

import tool_fig_config


parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input-dir', type=str, help='Input file', required=True)
parser.add_argument('--output', type=str, help='Output file', default="")
parser.add_argument('--varnames', type=str, nargs="+", help='Output file', default=["dMLTdt", "MLG_frc", "MLG_nonfrc"])
parser.add_argument('--time', type=str, nargs="+", help='Output file')
parser.add_argument('--thumbnail-offset', type=int, help='Output file', default=0)
parser.add_argument('--add-thumbnail-title', action="store_true")
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--no-sig', action="store_true")
parser.add_argument('--mark-box', type=str, nargs="*", default=[], choices=["strong_marine_warming"])

args = parser.parse_args()
print(args)

ds_stat = {}
for k in ["clim", "AR",]:
    ds_stat[k] = xr.open_dataset("%s/stat_%s.nc" % (args.input_dir, k))



plot_infos_scnario = {

    "clim" : {
        "title" : "All",
    },

    "AR" : {
        "title" : "AR",
    },

    "ARf" : {
        "title" : "AR free",
    },

    "AR-ARf" : {
        "title" : "AR minus AR free",
    }

}


G_scale = 0.5
G_scale2 = 0.2
G_scale3 = 0.1
G_scale4 = 0.01
plot_infos = {
    
    "dMLTdt" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{loc}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    },

    "MLG_frc" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{sfc}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_nonfrc" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{ocn}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_vdiff" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{vdiff}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_vdiff_entwep" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{ent}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_hdiff" : {
        "levels": np.linspace(-1, 1, 11) * G_scale3,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{hdif}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_vmix" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{vmix}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 


    "MLG_ent_wen" : {
        "levels": np.linspace(-1, 1, 11) * G_scale3,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{det}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 


    "MLG_ent" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{ent}} + \\dot{\\overline{\\Theta}}_{\\mathrm{det}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_adv" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\mathrm{adv}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

 
    "MLG_frc_sw" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{sw}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_frc_lw" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{lw}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_frc_sh" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{sen}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_frc_lh" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{lat}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_frc_fwf" : {
        "levels": np.linspace(-1, 1, 11) * G_scale2,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{fwf}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLG_frc_dilu" : {
        "levels": np.linspace(-1, 1, 11) * G_scale3,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "$ \\dot{\\overline{\\Theta}}_{\\mathrm{dilu}} $",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 


    "MLHADVT_g" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "Geo HADV",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "MLHADVT_ag" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "Ageo HADV",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 

    "ENT_ADV" : {
        "levels": np.linspace(-1, 1, 11) * G_scale,
        "levels_std": np.linspace(0, 2, 11),
        "label" : "Ent ADV",
        "unit"  : "$ 10^{-6} \\, \\mathrm{K} / \\mathrm{s} $",
        "factor" : 1e-6,
    }, 



    "MXLDEPTH" : {
        "levels": np.linspace(-1, 1, 11) * 20,
        "levels_std": np.linspace(0, 20, 21),
        "label" : "MXLDEPTH",
        "unit"  : "$\\mathrm{m}$",
        "factor" : 1.0,
    }, 
 
    "MLD" : {
        "levels": np.linspace(-1, 1, 11) * 20,
        "levels_std": np.linspace(0, 100, 21),
        "label" : "MLD",
        "unit"  : "$\\mathrm{m}$",
        "factor" : 1.0,
    }, 


    "IWV" : {
        "levels": np.linspace(0, 1, 11) * 50,
        "levels_std": np.linspace(0, 5, 11),
        "label" : "IWV",
        "unit"  : "$\\mathrm{kg} \\, / \\, \\mathrm{m}^2$",
        "factor" : 1.0,
    }, 

    "IVT" : {
        "levels": np.linspace(0, 1, 13) * 600,
        "levels_std": np.linspace(0, 1, 11) * 200,
        "label" : "IVT",
        "unit"  : "$\\mathrm{kg} \\, / \\, \\mathrm{m} \\, / \\, \\mathrm{s}$",
        "factor" : 1.0,
    }, 

    "SFCWIND" : {
        "levels": np.linspace(-1, 1, 13) * 6,
        "levels_std": np.linspace(0, 5, 11),
        "label" : "$\\left| \\vec{V}_{10} \\right|$",
        "unit"  : "$\\mathrm{m} \\, / \\, \\mathrm{s}$",
        "factor" : 1.0,
    }, 

    "lcc" : {
        "levels": np.linspace(-1, 1, 11) * 0.5,
        "levels_std": np.linspace(0, 1, 11) * 0.5,
        "label" : "LCC",
        "unit"  : "None",
        "factor" : 1.0,
    }, 

    "mcc" : {
        "levels": np.linspace(-1, 1, 11) * 0.5,
        "levels_std": np.linspace(0, 1, 11) * 0.5,
        "label" : "MCC",
        "unit"  : "None",
        "factor" : 1.0,
    }, 

    "hcc" : {
        "levels": np.linspace(-1, 1, 11) * 0.5,
        "levels_std": np.linspace(0, 1, 11) * 0.5,
        "label" : "HCC",
        "unit"  : "None",
        "factor" : 1.0,
    }, 


    "dTdz_b" : {
        "levels": np.linspace(-1, 1, 11) * 2,
        "levels_std": np.linspace(0, 1, 11) * 2,
        "label" : "$\\partial \\Theta_{\\eta - h} / \\partial z $",
        "unit"  : "$10^{-2} \\, \\mathrm{K} \\, / \\, \\mathrm{m}$",
        "factor" : 1e-2,
    }, 


}


# Plot data
print("Loading Matplotlib...")
import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    mpl.rc('font', size=15)
    mpl.rc('axes', labelsize=15)
     
 
  
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


print("done")

from scipy.stats import ttest_ind_from_stats

def student_t_test(mean1, std1, nobs1, mean2, std2, nobs2):

    shp = mean1.shape
    
    pvals = np.zeros(shp, dtype=np.float64)
    pvals[:] = np.nan
    
    tstats = pvals.copy()

    for j in range(shp[0]):
        for i in range(shp[1]):
            
            m1 = mean1[j, i]
            m2 = mean2[j, i]
            
            s1 = std1[j, i]
            s2 = std2[j, i]

            n1 = nobs1[j, i]
            n2 = nobs2[j, i]

            if np.all(np.isfinite([m1, m2, s1, s2, n1, n2])) and n1 >= 30 and n2 >= 30:
                
                result = ttest_ind_from_stats(
                    m1, s1, n1,
                    m2, s2, n2,
                    equal_var=False,
                    alternative='two-sided',
                )

                tstats[j, i] = result[0]
                pvals[j, i] = result[1]



    return tstats, pvals

cent_lon = 180.0

plot_lon_l = 100.0
plot_lon_r = 260.0
plot_lat_b = 10.0
plot_lat_t = 60.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

varnames = args.varnames

ncol = len(args.time)
nrow = len(varnames)

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 4.8,
    h = 2.0,
    wspace = 1.0,
    hspace = 0.5,
    w_left = 1.0,
    w_right = 2.2,
    h_bottom = 1.0,
    h_top = 1.0,
    ncol = ncol,
    nrow = nrow,
)


fig, ax = plt.subplots(
    nrow, ncol,
    figsize=figsize,
    subplot_kw=dict(projection=proj, aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
)

coords = ds_stat["clim"].coords
cmap = cm.get_cmap("bwr")

cmap.set_over("green")
cmap.set_under("yellow")

mappables = [ None for i in range(len(varnames)) ]

thumbnail_cnt = 0
for i, mon in enumerate(args.time):

    _ax = ax[:, i]
  
    # Set title for different month except for the whole average mon==7
    _ax[0].set_title(mon, size=30)

    for j, varname in enumerate(varnames):


        if varname == "BLANK":
            
            fig.delaxes(_ax[j])
            continue

        AR_da = ds_stat["AR"][varname].sel(time=mon)
        clim_da = ds_stat["clim"][varname].sel(time=mon)

        
        # 4 = annual_mean = the mean of annual means 
        # 8 = sampling count = the count of years whose annual means are used
        _mean1 = AR_da.sel(stat="mean").to_numpy()
        _mean2 = clim_da.sel(stat="mean").to_numpy()

        _std1 = AR_da.sel(stat="std").to_numpy()
        _std2 = clim_da.sel(stat="std").to_numpy()

        _nobs1 = AR_da.sel(stat="cnt").to_numpy()
        _nobs2 = clim_da.sel(stat="cnt").to_numpy()

        _, pvals = student_t_test(_mean1, _std1, _nobs1, _mean2, _std2, _nobs2)        
        

        plot_info = plot_infos[varname]

        if varname in ["IWV", "IVT"]:

            _plot_data = ( AR_da.sel(stat="mean").to_numpy() + clim_da.sel(stat="mean").to_numpy() ) / plot_info["factor"]
            mappables[j] = _ax[j].contourf(
                coords["lon"], coords["lat"],
                _plot_data,
                levels=plot_info["levels"],
                cmap="GnBu",
                extend="max",
                transform=proj_norm,
            )

        else:
            _plot_data = AR_da.sel(stat="mean").to_numpy() / plot_info["factor"]
            mappables[j] = _ax[j].contourf(
                coords["lon"], coords["lat"],
                _plot_data,
                levels=plot_info["levels"],
                cmap=cmap, 
                extend="both", 
                transform=proj_norm,
            )
       
        if args.add_thumbnail_title :
            _ax[j].set_title("(%s)" % ("abcdefghijklmnopqrstuvwxyz"[thumbnail_cnt + args.thumbnail_offset]))
            thumbnail_cnt += 1

        
        if varname == "SFCWIND":
            
            readjust_U = (ds_stat["AR"]["u10"].sel(time=mon, stat="mean") / np.cos(np.deg2rad(coords["lat"]))).to_numpy()
            readjust_V =  ds_stat["AR"]["v10"].sel(time=mon, stat="mean").to_numpy()
            _ax[j].streamplot(coords["lon"], coords["lat"], readjust_U, readjust_V, transform=proj_norm, density=1, color="dodgerblue")
            
        if not args.no_sig:
            # Plot the hatch to denote significant data
            _dot = np.zeros_like(_plot_data)
            _significant_idx =  (pvals <= 0.05) 

            _dot[ _significant_idx                 ] = 0.75
            _dot[ np.logical_not(_significant_idx) ] = 0.25

            # Remove insignificant data
            #_diff[np.logical_not(_significant_idx)] = np.nan

            cs = _ax[j].contourf(coords["lon"], coords["lat"], _dot, colors='none', levels=[0, 0.5, 1], hatches=[None, "."], transform=proj_norm)

            # Remove the contour lines for hatches 
            for _, collection in enumerate(cs.collections):
                collection.set_edgecolor((.2, .2, .2))
                collection.set_linewidth(0.)

        # Plot the standard deviation
        std = AR_da.sel(stat="std") / plot_info["factor"]
        cs = _ax[j].contour(coords["lon"], coords["lat"], std, levels=plot_info["levels_std"], colors="k", linestyles='-',linewidths=1, transform=proj_norm, alpha=0.8, zorder=10)


        fmt = "%d" if np.all( np.abs( np.array(plot_info["levels_std"]) % 1 ) < 1e-5 ) else "%.1f"
        _ax[j].clabel(cs, fmt=fmt)


        if ("strong_marine_warming" in args.mark_box) and (i == 0) and (j == 0):
            lon_rng = np.array([-135, -120])  % 360
            lat_rng = np.array([15, 25])
            _ax[j].add_patch(
                Rectangle(
                    (lon_rng[0], lat_rng[0]), lon_rng[1] - lon_rng[0], lat_rng[1] - lat_rng[0],
                    linewidth=2, edgecolor="cyan", facecolor='none', transform=proj_norm,
                )
            )
            
    for __ax in _ax: 

        __ax.set_global()
        #__ax.gridlines()
        __ax.coastlines(color='gray')
        __ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

        gl = __ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top   = False
        gl.ylabels_right = False

        #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
        gl.xlocator = mticker.FixedLocator([120, 150, 180, -150, -120])#np.arange(-180, 181, 30))
        gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])
        
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12, 'color': 'black'}
        gl.ylabel_style = {'size': 12, 'color': 'black'}



for i, mappable in enumerate(mappables):
   
    if mappable is None:
        continue 

    cax = tool_fig_config.addAxesNextToAxes(fig, ax[i, -1], "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.00)

    plot_info = plot_infos[varnames[i]]
    unit_str = "" if plot_info["unit"] == "" else " [ %s ]" % (plot_info["unit"],)
    cb.ax.set_ylabel("%s\n%s" % (plot_info["label"], unit_str), size=25)


if not args.no_display:
    plt.show()

if args.output != "":

    print("Saving output: ", args.output)    
    fig.savefig(args.output, dpi=200)


print("Finished.")


