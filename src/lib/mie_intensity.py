from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from lib.mie_scattering import mie_scattering

if TYPE_CHECKING:
    from numpy import ndarray


@dataclass
class MieIntensityResults:
    wl: float
    mat: str
    r: float
    n_med: float
    theta: "ndarray[float]"
    i_tm: "ndarray[float]"
    i_te: "ndarray[float]"


def i_parallel(a: "ndarray", b: "ndarray", pi: "ndarray", tau: "ndarray") -> float:
    """
    Calculate the parallel intensity component.

    Parameters
    ----------
    a : ndarray
        Mie coefficients a_n.
    b : ndarray
        Mie coefficients b_n.
    pi : ndarray
        Angular function pi_n.
    tau : ndarray
        Angular function tau_n.

    Returns
    -------
    float
        Parallel intensity component.
    """
    max_n = min(len(a), len(b), len(pi), len(tau))

    return (
        np.absolute(
            sum(
                (2 * n + 1)
                / (n * (n + 1))
                * (a[n - 1] * tau[n - 1] + b[n - 1] * pi[n - 1])
                for n in range(1, max_n + 1)
            )
        )
        ** 2
    )


def i_perpendicular(a: "ndarray", b: "ndarray", pi: "ndarray", tau: "ndarray") -> float:
    """
    Calculate the perpendicular intensity component.

    Parameters
    ----------
    a : ndarray
        Mie coefficients a_n.
    b : ndarray
        Mie coefficients b_n.
    pi : ndarray
        Angular function pi_n.
    tau : ndarray
        Angular function tau_n.

    Returns
    -------
    float
        Perpendicular intensity component.
    """
    max_n = min(len(a), len(b), len(pi), len(tau))

    return (
        np.absolute(
            sum(
                (2 * n + 1)
                / (n * (n + 1))
                * (a[n - 1] * pi[n - 1] + b[n - 1] * tau[n - 1])
                for n in range(1, max_n + 1)
            )
        )
        ** 2
    )


def mie_intensity(wl: float, mat: str, r: float, n_med: float) -> MieIntensityResults:
    """
    Calculates Mie intensity for a spherical nanoparticle.

    Parameters
    ----------
    wl : float
        Incident wavelength in nanometers.
    mat : str
        Material symbol (e.g. 'Si', 'Au', 'Ag', etc.).
    r : float
        Nanoparticle radius in nanometers.
    n_med : float
        Refractive index of the surrounding medium.

    Returns
    -------
    MieIntensityResults
        Contains wavelength, material, radius, medium refractive index, angles, and intensities.
    """
    # === Constants ===
    theta = np.arange(0, 2 * np.pi, 0.01)

    # Angular functions
    pi1 = 1
    pi2 = 3 * np.cos(theta)
    tau1 = np.cos(theta)
    tau2 = 6 * np.cos(theta) ** 2 - 3

    # === Mie coefficients a1, b1, a2, b2 ===
    mie_result = mie_scattering(
        wl_range=np.array([wl]),
        mat=mat,
        r=r,
        n_med=n_med,
    )

    a1 = mie_result.a1[0]
    b1 = mie_result.b1[0]
    a2 = mie_result.a2[0]
    b2 = mie_result.b2[0]

    # === Scattering intensities ===
    I_TM = i_parallel([a1, a2], [b1, b2], [pi1, pi2], [tau1, tau2])
    I_TE = i_perpendicular([a1, a2], [b1, b2], [pi1, pi2], [tau1, tau2])

    return MieIntensityResults(
        wl=wl,
        mat=mat,
        r=r,
        n_med=n_med,
        theta=theta,
        i_tm=I_TM,
        i_te=I_TE,
    )
