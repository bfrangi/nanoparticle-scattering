# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 10:28:18 2025

@author: BRGARCIA
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, spherical_yn

def Mie_theory_intensity(mater):
    """
    Python translation of:
    Mie_theory_Calculations_Intensity_MiPhot.m

    Parameters
    ----------
    mater : str
        Symbol of the material ('Si', 'Ge', 'Au', ...)

    Output:
        Saves ScatterIntensity.dat with [theta, I_TM, I_TE]
        Produces a polar plot of scattering intensities
    """

    # -----------------------------
    # User input
    # -----------------------------
    r = float(input("Enter the value of the radius of the nanoparticle in nm: "))
    lambda_in = float(input("Enter the incident wavelength in nm: "))
    nmed = float(input("Enter the value of the refractive index of the surrounding medium: "))

    # -----------------------------
    # Constants
    # -----------------------------
    w = 1
    theta = np.arange(0, 2*np.pi, 0.01)

    # Angular functions
    pi1 = 1
    pi2 = 3*np.cos(theta)
    tau1 = np.cos(theta)
    tau2 = 6*np.cos(theta)**2 - 3

    # -----------------------------
    # Load optical constants
    # -----------------------------
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

    mu = 1

    # -----------------------------
    # Linear interpolation for optical constants
    # -----------------------------
    ni = np.interp(lambda_in, lambdaMat[::-1], nMat[::-1])
    ki = np.interp(lambda_in, lambdaMat[::-1], kiMat[::-1])

    # -----------------------------
    # Mie calculations
    # -----------------------------
    k = 2*np.pi*nmed / lambda_in      # wavenumber
    x = k * r
    m = complex(ni, ki) / nmed
    m1 = m / mu
    m1b = 1/m1
    z = m * x

    # spherical Bessel and Hankel functions
    def sph_bessel_j(n, val):
        return np.sqrt(np.pi/(2*val)) * spherical_jn(n, val)

    def sph_bessel_y(n, val):
        return np.sqrt(np.pi/(2*val)) * spherical_yn(n, val)

    # j_n(x)
    bj0 = sph_bessel_j(0.5, x)
    bj1 = sph_bessel_j(1.5, x)
    bj2 = sph_bessel_j(2.5, x)

    # y_n(x)
    by0 = sph_bessel_y(0.5, x)
    by1 = sph_bessel_y(1.5, x)
    by2 = sph_bessel_y(2.5, x)

    # h_n(x)
    bh0 = bj0 + 1j*by0
    bh1 = bj1 + 1j*by1
    bh2 = bj2 + 1j*by2

    # derivatives
    d_bj  = x*bj0 - w*bj1
    d_bh  = x*bh0 - w*bh1
    d_bj2 = x*bj1 - 2*bj2
    d_bh2 = x*bh1 - 2*bh2

    # for z
    bj0z = sph_bessel_j(0.5, z)
    bj1z = sph_bessel_j(1.5, z)
    bj2z = sph_bessel_j(2.5, z)

    by0z = sph_bessel_y(0.5, z)
    by1z = sph_bessel_y(1.5, z)
    by2z = sph_bessel_y(2.5, z)

    bh0z = bj0z + 1j*by0z
    bh1z = bj1z + 1j*by1z
    bh2z = bj2z + 1j*by2z

    d_bjz  = z*bj0z - w*bj1z
    d_bjz2 = z*bj1z - 2*bj2z

    # -----------------------------
    # Mie coefficients a1, b1, a2, b2
    # -----------------------------
    a1 = (m1*(z*bj1z*d_bj) - d_bjz*(x*bj1)) / (m1*(z*bj1z*d_bh) - d_bjz*(x*bh1))
    b1 = (m1b*(z*bj1z*d_bj) - d_bjz*(x*bj1)) / (m1b*(z*bj1z*d_bh) - d_bjz*(x*bh1))

    a2 = (m1*(z*bj2z*d_bj2) - d_bjz2*(x*bj2)) / (m1*(z*bj2z*d_bh2) - d_bjz2*(x*bh2))
    b2 = (m1b*(z*bj2z*d_bj2) - d_bjz2*(x*bj2)) / (m1b*(z*bj2z*d_bh2) - d_bjz2*(x*bh2))

    # -----------------------------
    # Scattering intensities
    # -----------------------------
    # Write the expressions of the parallel and a perpendicular polarizaed intensities
    I_TM = 
    I_TE = 

    # -----------------------------
    # Output plot
    # -----------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.plot(theta, I_TM, label='TM')
    ax.plot(theta, I_TE, label='TE')
    ax.legend()
    plt.show()

    # -----------------------------
    # Save data
    # -----------------------------
    out = np.column_stack([theta, I_TM, I_TE])
    np.savetxt("ScatterIntensity.dat", out)

    print("\nScatterIntensity.dat saved successfully.\n")
