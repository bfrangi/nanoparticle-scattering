import os
from pathlib import Path


def get_root_path() -> str:
    """
    Get the root project path.

    Returns
    -------
    str
        The root project path.
    """
    script_dir = os.path.dirname(__file__)
    return Path(os.path.normpath(os.path.join(script_dir, "../..")))


def get_figures_path() -> str:
    """
    Get the figures folder path.

    Returns
    -------
    str
        The figures folder path.
    """
    return get_root_path() / "figures"
