from typing import TYPE_CHECKING

import numpy as np
from matplotlib import pyplot as plt

from lib.files import get_figures_path
from lib.mie_scattering import mie_scattering
from lib.plots import article_dpi, article_figsize, article_tight, get_cmap, use_latex

if TYPE_CHECKING:
    from lib.mie_scattering import MieScatteringResults

materials = [
    "Au",
    "Ag",
]
n_med = 1.333
r = 30  # in nm
min_wl = 300  # in nm
max_wl = 900  # in nm
steps = 1000

wl_range = np.linspace(min_wl, max_wl, steps)

results: "list[MieScatteringResults]" = []
cmap = get_cmap("viridis", 3)
use_latex()

for mat in materials:
    mie_results = mie_scattering(
        wl_range=wl_range,
        mat=mat,
        r=r,
        n_med=n_med,
    )
    results.append(mie_results)

for result in results:
    path = get_figures_path() / f"metal_extinction_{result.mat}.pdf"
    plt.figure(figsize=article_figsize)

    plt.plot(result.wl, result.Qext, label="Extinction", color=cmap(0))
    plt.plot(result.wl, result.Qsca, label="Scattering", color=cmap(1))
    plt.plot(result.wl, result.Qabs, label="Absorption", color=cmap(2))

    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Efficiency [a.u.]")
    legend = plt.legend(frameon=True, prop={"size": 7})
    legend.set_zorder(1001)
    legend.get_frame().set_edgecolor("none")
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.7)
    plt.tight_layout(**article_tight)
    plt.savefig(path, dpi=article_dpi)
    plt.clf()
