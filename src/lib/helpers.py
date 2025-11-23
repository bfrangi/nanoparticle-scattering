from typing import TYPE_CHECKING

from numpy import argmax

if TYPE_CHECKING:
    from numpy import ndarray


def get_x_of_max_y(x: "ndarray[float]", y: "ndarray[float]") -> float:
    """
    Get the x value at which the y is maximum.

    Parameters
    ----------
    x : ndarray[float]
        The wavelength data.
    y : ndarray[float]
        The transmission data.

    Returns
    -------
    float
        The x value at which y is maximum.
    """
    return x[argmax(y)]
