import matplotlib.pyplot as plt

article_tight = {"pad": 0.1, "rect": (0.01, 0.01, 0.99, 0.99)}
article_scale = 0.3
article_figsize = (8 * article_scale, 6 * article_scale)
article_dpi = 300

cmaps = {
    "brg": 2,
    "CMRmap": 2,
    "coolwarm": 1,
    "copper": 1,
    "twilight": 4 / 3,
    "viridis": 3 / 2,
    "autumn": 3 / 2,
    "berlin": 1,
    "gnuplot": 4 / 3,
    "summer": 3 / 2,
    "winter": 1,
}


def get_cmap(name: str, n: int):
    """Get a colormap with n distinct colors.

    Parameters
    ----------
    name : str
        The name of the colormap.
    n : int
        The number of distinct colors needed.

    Returns
    -------
    matplotlib.colors.Colormap
        The requested colormap with n distinct colors.
    """
    if name in cmaps:
        factor = cmaps[name]
        n = int(n * factor)
    return plt.get_cmap(name, n)


def use_latex():
    """Use LaTeX for rendering text in plots."""
    plt.rcParams.update({"text.usetex": True, "font.family": "Computer Modern"})
