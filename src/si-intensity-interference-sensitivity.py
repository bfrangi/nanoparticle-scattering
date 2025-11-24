from typing import TYPE_CHECKING

import numpy as np
from matplotlib import pyplot as plt

from lib.files import get_figures_path
from lib.mie_intensity import mie_intensity
from lib.plots import (
    article_dpi,
    article_figsize_polar,
    article_tight,
    get_cmap,
    use_latex,
)

if TYPE_CHECKING:
    from lib.mie_intensity import MieIntensityResults

mat = "Si"
n_med = 1.333
r = 100  # in nm
wls = [
    423.4932349323493,
    430.29730297302973,
    467.119671196712,
    504.8360483604836,
    562.5206252062521,
    725.6502565025651,
    840.6174061740618,
]  # in nm

results: "list[list[MieIntensityResults]]" = []
use_latex()

for wl in wls:
    results_this = []
    wls_this = wl * np.array([0.99, 1.0, 1.01])

    for wl_this in wls_this:
        result = mie_intensity(wl=wl_this, mat=mat, r=r, n_med=n_med)
        results_this.append(result)
    
    results.append(results_this)

for results_this, wl_center in zip(results, wls):
    cmap = get_cmap("viridis", len(results_this))
    path = (
        get_figures_path()
        / f"intensity_interference_sensitivity_{wl_center:.2f}_nm_{mat}.pdf"
    )

    fig = plt.figure(figsize=article_figsize_polar)
    ax = fig.add_subplot(111, polar=True)

    for i, result in enumerate(results_this):
        ax.plot(result.theta, result.i_tm, label=f"TM({result.wl:.2f} nm)", linestyle="-", color=cmap(i))
        ax.plot(result.theta, result.i_te, label=f"TE({result.wl:.2f} nm)", linestyle="--", color=cmap(i))

    legend = plt.legend(loc=(-0.15, -0.15), frameon=True, prop={"size": 7})
    legend.set_zorder(1001)
    legend.get_frame().set_edgecolor("none")
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.7)
    plt.tight_layout(**article_tight)
    plt.savefig(path, dpi=article_dpi)
    plt.clf()
