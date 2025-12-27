from typing import TYPE_CHECKING

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
    564.5606456064561,
    839.789397893979,
]  # in nm

results: "list[MieIntensityResults]" = []
cmap = get_cmap("viridis", 2)
use_latex()

for wl in wls:
    result = mie_intensity(wl=wl, mat=mat, r=r, n_med=n_med)
    results.append(result)

for result in results:
    path = (
        get_figures_path()
        / f"intensity_interference_{result.wl:.2f}_nm_{result.mat}.pdf"
    )

    fig = plt.figure(figsize=article_figsize_polar)
    ax = fig.add_subplot(111, polar=True)
    ax.plot(result.theta, result.i_tm, label="TM", color=cmap(0))
    ax.plot(result.theta, result.i_te, label="TE", color=cmap(1))
    legend = plt.legend(frameon=True, prop={"size": 7})
    legend.set_zorder(1001)
    legend.get_frame().set_edgecolor("none")
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.7)
    plt.tight_layout(**article_tight)
    plt.savefig(path, dpi=article_dpi)
    plt.clf()
