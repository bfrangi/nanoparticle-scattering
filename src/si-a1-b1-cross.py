import numpy as np

from lib.mie_scattering import mie_scattering

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

wls = []

for i, (a1, b1, wl) in enumerate(zip(mie_result.a1, mie_result.b1, mie_result.wl)):
    a1_i = np.real(a1)
    b1_i = np.real(b1)
    a1_prev = np.real(mie_result.a1[i - 1]) if i > 0 else a1_i
    b1_prev = np.real(mie_result.b1[i - 1]) if i > 0 else b1_i

    if (a1_i - b1_i) * (a1_prev - b1_prev) < 0:
        wls.append(wl)

print("Wavelengths where a1 crosses b1:")
for wl in wls:
    print(f"{wl} nm")
