import numpy as np
from matplotlib import pyplot as plt

from lib.files import get_figures_path
from lib.mie_scattering import mie_scattering
from lib.plots import (
    article_dpi,
    article_figsize_wide,
    article_tight,
    get_cmap,
    use_latex,
)

mat = "Si"
n_med = 1.333
r = 100  # in nm
min_wl = 300  # in nm
max_wl = 900  # in nm
steps = 1000

wl_range = np.linspace(min_wl, max_wl, steps)

cmap = get_cmap("viridis", 5)
use_latex()

mie_result = mie_scattering(
    wl_range=wl_range,
    mat=mat,
    r=r,
    n_med=n_med,
)

path = get_figures_path() / f"extinction_breakdown_{mie_result.mat}.pdf"
plt.figure(figsize=article_figsize_wide)

plt.plot(mie_result.wl, mie_result.Qext, label="Extinction", color=cmap(0))

plt.plot(
    mie_result.wl,
    np.real(mie_result.a1),
    label="Re$(a_1)$",
    linestyle="--",
    color=cmap(1),
    zorder=10,
)
plt.annotate(
    "ed",
    xy=(620, 1),
    xytext=(630, 2),
    arrowprops=dict(arrowstyle="->", color=cmap(1), shrinkA=0, shrinkB=0),
    color=cmap(1),
    zorder=10,
)
plt.annotate(
    "ed",
    xy=(425, 0.7),
    xytext=(425, 1.4),
    arrowprops=dict(arrowstyle="->", color=cmap(1), shrinkA=0, shrinkB=0),
    color=cmap(1),
    zorder=10,
)

plt.plot(
    mie_result.wl,
    np.real(mie_result.a2),
    label="Re$(a_2)$",
    linestyle="--",
    color=cmap(2),
    zorder=20,
)
plt.annotate(
    "eq",
    xy=(505, 0.74),
    xytext=(530, 1),
    arrowprops=dict(arrowstyle="->", color=cmap(2), shrinkA=0, shrinkB=0),
    color=cmap(2),
    zorder=20,
)
plt.annotate(
    "eq",
    xy=(405, 0.58),
    xytext=(410, 1.9),
    arrowprops=dict(arrowstyle="->", color=cmap(2), shrinkA=0, shrinkB=0),
    color=cmap(2),
    zorder=20,
)

plt.plot(
    mie_result.wl,
    np.real(mie_result.b1),
    label="Re$(b_1)$",
    linestyle="--",
    color=cmap(3),
    zorder=30,
)
plt.annotate(
    "md",
    xy=(775, 1.02),
    xytext=(765, 1.5),
    arrowprops=dict(arrowstyle="->", color=cmap(3), shrinkA=0, shrinkB=0),
    color=cmap(3),
    zorder=30,
)
plt.annotate(
    "md",
    xy=(455, 0.78),
    xytext=(500, 1.4),
    arrowprops=dict(arrowstyle="->", color=cmap(3), shrinkA=0, shrinkB=0),
    color=cmap(3),
    zorder=30,
)

plt.plot(
    mie_result.wl,
    np.real(mie_result.b2),
    label="Re$(b_2)$",
    linestyle="--",
    color=cmap(4),
    zorder=40,
)
plt.annotate(
    "mq",
    xy=(580, 0.78),
    xytext=(580, 1.2),
    arrowprops=dict(arrowstyle="->", color=cmap(4), shrinkA=0, shrinkB=0),
    color=cmap(4),
    zorder=40,
)
plt.annotate(
    "mq",
    xy=(420, 0.36),
    xytext=(480, 1.9),
    arrowprops=dict(arrowstyle="->", color=cmap(4), shrinkA=0, shrinkB=0),
    color=cmap(4),
    zorder=40,
)


plt.xlabel("Wavelength [nm]")
plt.ylabel("Eff. \& Mie coeff. [a.u.]")
legend = plt.legend(frameon=True, prop={"size": 7})
legend.set_zorder(1001)
legend.get_frame().set_edgecolor("none")
legend.get_frame().set_facecolor("white")
legend.get_frame().set_alpha(0.7)
plt.tight_layout(**article_tight)
plt.savefig(path, dpi=article_dpi)
plt.clf()
