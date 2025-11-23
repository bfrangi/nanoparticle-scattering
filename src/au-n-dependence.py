from typing import TYPE_CHECKING

import numpy as np
from matplotlib import pyplot as plt

from lib.files import get_figures_path
from lib.mie_scattering import mie_scattering
from lib.plots import article_dpi, article_figsize, article_tight, get_cmap, use_latex
from lib.helpers import get_x_of_max_y

if TYPE_CHECKING:
    from lib.mie_scattering import MieScatteringResults

material = "Au"
ns_med = [1, 1.2, 1.4, 1.6, 1.8, 2]
r = 30  # in nm
min_wl = 350  # in nm
max_wl = 750  # in nm
steps = 1000

wl_range = np.linspace(min_wl, max_wl, steps)

results: "list[MieScatteringResults]" = []
cmap = get_cmap("viridis", len(ns_med))
use_latex()

for n_med in ns_med:
    mie_results = mie_scattering(
        wl_range=wl_range,
        mat=material,
        r=r,
        n_med=n_med,
    )
    results.append(mie_results)

path = get_figures_path() / "n_med_dependence_extinction.pdf"
plt.figure(figsize=article_figsize)
for i, result in enumerate(results):
    plt.plot(result.wl, result.Qext, label=f"$n = {result.n_med}$", color=cmap(i))
plt.xlabel("Wavelength [nm]")
ylab = plt.ylabel("Extinction efficiency [a.u.]")
ylab.set_position((ylab.get_position()[0], 0.4))
legend = plt.legend(frameon=True, prop={"size": 7})
legend.set_zorder(1001)
legend.get_frame().set_edgecolor("none")
legend.get_frame().set_facecolor("white")
legend.get_frame().set_alpha(0.7)
plt.tight_layout(**article_tight)
plt.savefig(path, dpi=article_dpi)
plt.clf()


peak_ext_eff = [result.Qext.max() for result in results]
path = get_figures_path() / "n_med_dependence_max_extinction.pdf"
plt.figure(figsize=article_figsize)
plt.plot(ns_med, peak_ext_eff, color=cmap(0))
plt.xlabel("Refractive index [a.u.]")
ylab = plt.ylabel("Peak extinction eff. [a.u.]")
ylab.set_position((ylab.get_position()[0], 0.4))
plt.tight_layout(**article_tight)
plt.savefig(path, dpi=article_dpi)
plt.clf()

peak_ext_pos = [get_x_of_max_y(result.wl, result.Qext) for result in results]
path = get_figures_path() / "n_med_dependence_pos_max_extinction.pdf"
plt.figure(figsize=article_figsize)
plt.plot(ns_med, peak_ext_pos, color=cmap(0))
plt.xlabel("Refractive index [a.u.]")
ylab = plt.ylabel("Peak extinction pos. [nm]")
ylab.set_position((ylab.get_position()[0], 0.4))
plt.tight_layout(**article_tight)
plt.savefig(path, dpi=article_dpi)
plt.clf()
