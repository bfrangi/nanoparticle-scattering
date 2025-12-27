import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema

from lib.files import get_figures_path
from lib.mie_scattering import mie_scattering
from lib.plots import (
    article_dpi,
    article_figsize_stacked,
    article_tight,
    get_cmap,
    use_latex,
)

mat = "Si"
n_med = 1.333
r = 100  # in nm
min_wl = 300  # in nm
max_wl = 900  # in nm
steps = 100000

wl_range = np.linspace(min_wl, max_wl, steps)

mie_result = mie_scattering(
    wl_range=wl_range,
    mat=mat,
    r=r,
    n_med=n_med,
)

cmap = get_cmap("viridis", 3)
use_latex()
fig, (ax1, ax2) = plt.subplots(
    2,
    1,
    sharex=True,
    figsize=article_figsize_stacked,
)

ax1.plot(
    mie_result.wl,
    np.real(mie_result.a1),
    c=cmap(0),
    linestyle="-",
    label="Re($a_1$)",
)
ax1.plot(
    mie_result.wl,
    np.imag(mie_result.a1),
    c=cmap(0),
    linestyle="--",
    label="Im($a_1$)",
)
ax1.plot(
    mie_result.wl,
    np.real(mie_result.b1),
    c=cmap(2),
    linestyle="-",
    label="Re($b_1$)",
)
ax1.plot(
    mie_result.wl,
    np.imag(mie_result.b1),
    c=cmap(2),
    linestyle="--",
    label="Im($b_1$)",
)
ax1.set_ylabel("Mie coeff. [a.u.]")

proximity = np.abs(mie_result.a1 - mie_result.b1)
threshold = 0.1

min_indices = argrelextrema(proximity, np.less)[0]
min_indices = min_indices[proximity[min_indices] < threshold]
min_values = proximity[min_indices]
min_wls = mie_result.wl[min_indices]

ax2.plot(mie_result.wl, proximity, c=cmap(1))
ax2.set_xlabel("Wavelength [nm]")
ax2.set_ylabel("$\mid a_1 - b_1\mid$ [a.u.]")
ax2.set_yscale("log")

for val, wl in zip(min_values, min_wls):
    print(f"{wl} nm")
    ax2.annotate(
        f"{wl:.2f} nm",
        xy=(wl, val),
        xytext=(wl - 250, val * 2),
        arrowprops=dict(arrowstyle="->", color=cmap(1), shrinkA=0, shrinkB=0),
        color=cmap(1),
        zorder=10,
    )

legend = ax1.legend(frameon=True, prop={"size": 7})
legend.set_zorder(1001)
legend.get_frame().set_edgecolor("none")
legend.get_frame().set_facecolor("white")
legend.get_frame().set_alpha(0.7)

plt.tight_layout(**article_tight)
path = get_figures_path() / "si_a1_b1_cross.pdf"
plt.savefig(path, dpi=article_dpi)