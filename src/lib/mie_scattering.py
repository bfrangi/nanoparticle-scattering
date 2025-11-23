from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.special import spherical_jn, spherical_yn

from lib.files import get_root_path

if TYPE_CHECKING:
    from numpy.typing import NDArray


@dataclass
class MieScatteringResults:
    wl: "NDArray[np.float64]"
    mat: str
    r: float
    n_med: float
    a1: "NDArray[np.complex128]"
    b1: "NDArray[np.complex128]"
    a2: "NDArray[np.complex128]"
    b2: "NDArray[np.complex128]"
    Qext: "NDArray[np.float64]"
    Qsca: "NDArray[np.float64]"
    Qabs: "NDArray[np.float64]"


def c_sca(k: float, a: "NDArray[np.complex128]", b: "NDArray[np.complex128]") -> float:
    """
    Calculate the scattering cross section.

    Parameters
    ----------
    k : float
        Wave number.
    a : np.ndarray
        Mie coefficients a_n.
    b : np.ndarray
        Mie coefficients b_n.

    Returns
    -------
    float
        Scattering cross section.
    """
    max_n = min(len(a), len(b))

    return (2 / k**2) * sum(
        (2 * n + 1) * (np.absolute(a[n - 1]) ** 2 + np.absolute(b[n - 1]) ** 2)
        for n in range(1, max_n + 1)
    )


def c_ext(k: float, a: "NDArray[np.complex128]", b: "NDArray[np.complex128]") -> float:
    """
    Calculate the extinction cross section.

    Parameters
    ----------
    k : float
        Wave number.
    a : np.ndarray
        Mie coefficients a_n.
    b : np.ndarray
        Mie coefficients b_n.

    Returns
    -------
    float
        Extinction cross section.
    """
    max_n = min(len(a), len(b))

    return (2 / k**2) * sum(
        (2 * n + 1) * np.real(a[n - 1] + b[n - 1]) for n in range(1, max_n + 1)
    )


def mie_scattering(
    wl_range: "NDArray[np.float64]", mat: str, r: float, n_med: float
) -> MieScatteringResults:
    """
    Calculates Mie coefficients and optical efficiencies for a spherical nanoparticle.

    Parameters
    ----------
    wl_range : NDArray[np.float64]
        Range of incident wavelengths in nanometers.
    mat : str
        Material symbol (e.g. 'Si', 'Au', 'Ag', etc.).
    r : float
        Nanoparticle radius in nanometers.
    n_med : float
        Refractive index of the surrounding medium.

    Returns
    -------
    MieScatteringResults
        Contains wavelength, material, radius, medium refractive index, Mie coefficients, and efficiencies.
    """

    print(f"\nRunning Mie theory calculations for {mat} nanoparticles...")

    # === Load material data ===
    # Material files must contain: wavelength(µm), n(real part), k(imag part)
    path = get_root_path() / "data/refractive-index"

    file_map = {
        "Si": "Si.txt",
        "Ge": "Ge.txt",
        "GaP": "GaP.txt",
        "GaAs": "GaAs.txt",
        "AlSb": "AlSb.txt",
        "AlAs": "AlAs.txt",
        "TiO2": "TiO2.txt",
        "Au": "Au.txt",
        "Ag": "Ag.txt",
    }

    if mat not in file_map:
        raise ValueError("Unknown material name. Check spelling.")

    filename = path / file_map[mat]
    Mat_OC = np.loadtxt(filename)

    # Convert µm → nm
    lambdaMat = Mat_OC[:, 0] * 1e3
    nMat = Mat_OC[:, 1]
    kiMat = Mat_OC[:, 2]

    mu = 1.0  # Magnetic permeability of vacuum

    # === Linear interpolation of optical constants ===
    ni = np.interp(wl_range, lambdaMat[::-1], nMat[::-1])
    ki = np.interp(wl_range, lambdaMat[::-1], kiMat[::-1])

    # === Initialize result arrays ===
    a1 = np.zeros_like(wl_range, dtype=complex)
    b1 = np.zeros_like(wl_range, dtype=complex)
    a2 = np.zeros_like(wl_range, dtype=complex)
    b2 = np.zeros_like(wl_range, dtype=complex)

    Qext = np.zeros_like(wl_range, dtype=float)
    Qsca = np.zeros_like(wl_range, dtype=float)
    Qabs = np.zeros_like(wl_range, dtype=float)

    w = 1  # Constant used in the derivative expressions

    # === Main Mie calculations ===
    for j, lam in enumerate(wl_range):
        k = 2 * np.pi * n_med / lam  # wave number
        x = k * r  # size parameter

        m = (ni[j] + 1j * ki[j]) / n_med  # complex refractive index
        m1 = m / mu
        m1b = 1 / m1
        z = m * x

        # Spherical Bessel functions (order 0, 1, 2)
        bj = [spherical_jn(n, x) for n in range(3)]
        by = [spherical_yn(n, x) for n in range(3)]
        bh = [bj[n] + 1j * by[n] for n in range(3)]

        bjz = [spherical_jn(n, z) for n in range(3)]

        # Derivatives
        d_bj = [x * bj[0] - w * bj[1], x * bj[1] - 2 * bj[2]]
        d_bh = [x * bh[0] - w * bh[1], x * bh[1] - 2 * bh[2]]
        d_bjz = [z * bjz[0] - w * bjz[1], z * bjz[1] - 2 * bjz[2]]

        # Mie coefficients
        a1[j] = (m1 * (z * bjz[1] * d_bj[0]) - (d_bjz[0] * x * bj[1])) / (
            m1 * (z * bjz[1] * d_bh[0]) - (d_bjz[0] * x * bh[1])
        )
        b1[j] = (m1b * (z * bjz[1] * d_bj[0]) - (d_bjz[0] * x * bj[1])) / (
            m1b * (z * bjz[1] * d_bh[0]) - (d_bjz[0] * x * bh[1])
        )
        a2[j] = (m1 * (z * bjz[2] * d_bj[1]) - (d_bjz[1] * x * bj[2])) / (
            m1 * (z * bjz[2] * d_bh[1]) - (d_bjz[1] * x * bh[2])
        )
        b2[j] = (m1b * (z * bjz[2] * d_bj[1]) - (d_bjz[1] * x * bj[2])) / (
            m1b * (z * bjz[2] * d_bh[1]) - (d_bjz[1] * x * bh[2])
        )

        # Cross sections
        Cext = c_ext(k, [a1[j], a2[j]], [b1[j], b2[j]])
        Csca = c_sca(k, [a1[j], a2[j]], [b1[j], b2[j]])
        Cabs = Cext - Csca

        # Efficiencies
        Qext[j] = Cext / (np.pi * r**2)
        Qsca[j] = Csca / (np.pi * r**2)
        Qabs[j] = Cabs / (np.pi * r**2)

    return MieScatteringResults(
        wl=wl_range,
        mat=mat,
        r=r,
        n_med=n_med,
        a1=a1,
        b1=b1,
        a2=a2,
        b2=b2,
        Qext=Qext,
        Qsca=Qsca,
        Qabs=Qabs,
    )
