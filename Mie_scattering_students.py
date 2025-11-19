import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, spherical_yn

def Mie_Scattering(lambda_range, mater, r, nmed):
    """
    Calculates Mie coefficients and optical efficiencies for a spherical nanoparticle.

    Parameters
    ----------
    lambda_range : array-like
        Range of incident wavelengths in nanometers.
    mater : str
        Material symbol (e.g. 'Si', 'Au', 'Ag', etc.).
    r : float
        Nanoparticle radius in nanometers.
    nmed : float
        Refractive index of the surrounding medium.

    Output
    ------
    Saves two data files:
      - 'MieCoeff.dat' : wavelength and Mie coefficients (a1, b1, a2, b2)
      - 'Efficiencies.dat' : wavelength and efficiencies (Qext, Qsca, Qabs)
    Also plots both results.
    """

    print(f"\nRunning Mie theory calculations for {mater} nanoparticles...")

    # === Load material data ===
    # Material files must contain: wavelength(µm), n(real part), k(imag part)
    path = "G:/My Drive/DOCENCIA/Nanoelectronica y Nanofotonica/Sessions/Session 19 - Mie simulation/Materials/Refractive Index/"
    
    file_map = {
        'Si': 'Si.txt', 'Ge': 'Ge.txt', 'GaP': 'GaP.txt', 'GaAs': 'GaAs.txt',
        'AlSb': 'AlSb.txt', 'AlAs': 'AlAs.txt', 'TiO2': 'TiO2.txt',
        'Au': 'Au.txt', 'Ag': 'Ag.txt'
    }

    if mater not in file_map:
        raise ValueError("Unknown material name. Check spelling.")

    filename = path + file_map[mater]
    Mat_OC = np.loadtxt(filename)

    # Convert µm → nm
    lambdaMat = Mat_OC[:, 0] * 1e3
    nMat = Mat_OC[:, 1]
    kiMat = Mat_OC[:, 2]

    mu = 1.0  # Magnetic permeability of vacuum

    # === Linear interpolation of optical constants ===
    ni = np.interp(lambda_range, lambdaMat[::-1], nMat[::-1])
    ki = np.interp(lambda_range, lambdaMat[::-1], kiMat[::-1])

    # === Initialize result arrays ===
    a1 = np.zeros_like(lambda_range, dtype=complex)
    b1 = np.zeros_like(lambda_range, dtype=complex)
    a2 = np.zeros_like(lambda_range, dtype=complex)
    b2 = np.zeros_like(lambda_range, dtype=complex)

    Qext = np.zeros_like(lambda_range, dtype=float)
    Qsca = np.zeros_like(lambda_range, dtype=float)
    Qabs = np.zeros_like(lambda_range, dtype=float)

    w = 1  # Constant used in the derivative expressions

    # === Main Mie calculations ===
    for j, lam in enumerate(lambda_range):
        k = 2 * np.pi * nmed / lam          # wave number
        x = k * r                           # size parameter

        m = (ni[j] + 1j * ki[j]) / nmed     # complex refractive index
        m1 = m / mu
        m1b = 1 / m1
        z = m * x

        # Spherical Bessel functions (order 0, 1, 2)
        bj = [spherical_jn(n, x) for n in range(3)]
        by = [spherical_yn(n, x) for n in range(3)]
        bh = [bj[n] + 1j * by[n] for n in range(3)]

        bjz = [spherical_jn(n, z) for n in range(3)]
        byz = [spherical_yn(n, z) for n in range(3)]
        # bhz = [bjz[n] + 1j * byz[n] for n in range(3)]

        # Derivatives (following the same logic as MATLAB)
        d_bj = [x * bj[0] - w * bj[1],
                x * bj[1] - 2 * bj[2]]
        d_bh = [x * bh[0] - w * bh[1],
                x * bh[1] - 2 * bh[2]]
        d_bjz = [z * bjz[0] - w * bjz[1],
                 z * bjz[1] - 2 * bjz[2]]

        # Mie coefficients
        a1[j] = (m1 * (z * bjz[1] * d_bj[0]) - (d_bjz[0] * x * bj[1])) / \
                (m1 * (z * bjz[1] * d_bh[0]) - (d_bjz[0] * x * bh[1]))
        b1[j] = (m1b * (z * bjz[1] * d_bj[0]) - (d_bjz[0] * x * bj[1])) / \
                (m1b * (z * bjz[1] * d_bh[0]) - (d_bjz[0] * x * bh[1]))
        a2[j] = (m1 * (z * bjz[2] * d_bj[1]) - (d_bjz[1] * x * bj[2])) / \
                (m1 * (z * bjz[2] * d_bh[1]) - (d_bjz[1] * x * bh[2]))
        b2[j] = (m1b * (z * bjz[2] * d_bj[1]) - (d_bjz[1] * x * bj[2])) / \
                (m1b * (z * bjz[2] * d_bh[1]) - (d_bjz[1] * x * bh[2]))

        # Cross sections
        # Write here the correct expressions of the extinction, absroption and scattering cross sections
        Cext =        
        Csca = 
        Cabs = 

        # Efficiencies
        # Write here the correct expressions of the extinction, absorption and scattering efficiencies
        Qext[j] = 
        Qsca[j] = 
        Qabs[j] = 

    # === Plot efficiencies ===
    plt.figure(figsize=(8, 5))
    plt.plot(lambda_range, Qext, 'b', label='Extinction')
    plt.plot(lambda_range, Qsca, 'r', label='Scattering')
    plt.plot(lambda_range, Qabs, 'g', label='Absorption')
    plt.xlabel('Incident wavelength (nm)')
    plt.ylabel('Efficiency')
    plt.legend()
    plt.grid(True)
    plt.title(f'Optical Efficiencies - {mater}')
    plt.show()

    # === Plot Mie coefficients ===
    plt.figure(figsize=(8, 5))
    plt.plot(lambda_range, np.real(a1), label='a1')
    plt.plot(lambda_range, np.real(b1), label='b1')
    plt.plot(lambda_range, np.real(a2), label='a2')
    plt.plot(lambda_range, np.real(b2), label='b2')
    plt.xlabel('Incident wavelength (nm)')
    plt.ylabel('Mie Coefficients (Real part)')
    plt.legend()
    plt.grid(True)
    plt.title(f'Mie Coefficients - {mater}')
    plt.show()

    # === Save results ===
    np.savetxt("MieCoeff.dat", np.column_stack([lambda_range, a1.real, b1.real, a2.real, b2.real]))
    np.savetxt("Efficiencies.dat", np.column_stack([lambda_range, Qext, Qsca, Qabs]))

    print("\nFiles saved: 'MieCoeff.dat' and 'Efficiencies.dat'")

