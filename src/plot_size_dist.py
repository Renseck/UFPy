# =============================================================================
# This file functions to analyse the output of the *SMALL* SALSA2.0 Standalone model, which DOES NOT include
# atmospheric chemistry and cloud formation. There is no GitHub repository or documentation for this distribution.
# Investigation of the source code (driver.f90) seems to indicate that number concentration are
# currently written to output in num.dat.
# =============================================================================

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mc
import pandas as pd
import cmocean as co

def plot_size_dist(
    rdry, num, rows = [0], populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    ):

    ## make sure that row_nr is a list-like object
    try:
        iter(rows)
    except Exception:
        rows = [rows]

    nrows = 1
    ncols = len(populations)
    
    fig, axes = plt.subplots(
        nrows = nrows,
        ncols = ncols,
        figsize = (6*ncols, 4*nrows),
        sharex = True,
        sharey = True,
    )
    plt.tight_layout()
    fig.text(-0.01, 0.5, "# particles cm$^{-3}$", va = "center", rotation = "vertical")
    fig.text(0.5, 0.04, "Diameter (m)", ha = "center")

    for pop,ax in zip(populations, axes.ravel()):
        bins = [col for col in rdry.columns if pop in col]

        for n in rows:

            r_row = rdry.iloc[n][bins]
            N_row = num.iloc[n][bins]

            ax.plot(r_row, N_row, label=n)
        ax.set_title(f"Population {pop}")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(left=xmin, right=xmax)
        ax.set_ylim(bottom=ymin, top=ymax)
        if pop == "b":
            ax.legend(title = "Time", bbox_to_anchor = (1.2, 1.02))
    plt.show()


def define_bin_boundaries(populations = ['1a', '2a', '2b']):

    return {
        pop : (
            np.logspace(np.log10(3e-9), np.log10(50e-9), 4) if pop[0]=='1' else
            np.concatenate([
                np.logspace(np.log10(50e-9), np.log10(700e-9), 5, endpoint=True)[:-1],
                np.logspace(np.log10(700e-9), np.log10(10000e-9),4, endpoint=True),
            ])
        )
        for pop in populations
    }


def plot_size_dist_evolution(
    num, populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    vmin = None, vmax = None,
    ):

    bin_boundaries = define_bin_boundaries()

    nrows = len(populations)
    ncols = 1
    
    fig, axes = plt.subplots(
        nrows = nrows,
        ncols = ncols,
        figsize = (12*ncols, 4*nrows),
        sharex = True,
        sharey = True,
    )

    ## computing common color bar bounds:
    vmin = num.values.min() if vmin is None else vmin
    vmax = num.values.max() if vmax is None else vmax

    for pop,ax in zip(populations, axes.ravel()):
        bins = [col for col in rdry.columns if pop in col]

        ## combining bin boundaries (omitting double values)
        rbounds = [bounds for k,bounds in sorted(bin_boundaries.items()) if pop in k]
        rbounds = np.concatenate([rb[:-1] for rb in rbounds[:-1]]+[rbounds[-1]])

        ## time bounds
        tbounds = np.arange(num.shape[0]+1)

        ## generating meshgrid
        t,r = np.meshgrid(tbounds,rbounds)
      
        norm = mc.LogNorm(vmin = vmin, vmax = vmax)
        cls = ax.pcolormesh(t, r, num[bins].T, norm = norm, cmap = co.cm.dense)

        ax.set_title(f"Population {pop}")
        ax.set_yscale('log')
        ax.set_xlim(left=xmin, right=xmax)
        ax.set_ylim(bottom=ymin, top=ymax)
        
    fig.colorbar(cls, ax=axes.ravel().tolist())
    plt.savefig('size_distribution_LES_box.png', bbox_inches = 'tight', pad_inches = 0)
    # plt.close()
    plt.show()

        
if __name__ == '__main__':
    num  = pd.read_csv('../../HAM_box_OpenIFS/data/num.dat',  sep=r"\s+")
    rdry = pd.read_csv('../../HAM_box_OpenIFS/data/rdry.dat', sep=r"\s+")
    plot_size_dist(rdry, num, rows=[1,200,1000, 2000, 4000, 7080], ymin=1)
    plot_size_dist_evolution(num, vmin=1)
