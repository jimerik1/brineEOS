#!/usr/bin/env python3
"""
compare_water_density.py

Plot water density vs. pressure (0.1–100 MPa) at several temperatures
using:
  • A simple polynomial+compressibility model
  • IAPWS-IF97 Region 1 via the iapws package
"""

import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97

def rho_simple(T, P):
    """
    Simple polynomial + compressibility model.
    T in K, P in MPa → ρ in kg/m³
    """
    T0, P0, rho0 = 273.15, 0.1, 999.975
    alpha, beta = 0.000214, 0.46e-9
    dT = T - T0
    dP = (P - P0) * 1e6
    return rho0 * (1 + beta*dP) / (1 + alpha*dT)

def rho_iapws(T, P):
    """
    IAPWS-IF97 Region 1 density via iapws.
    Clamps P to [0.1, 100] MPa.
    """
    P_clamped = np.clip(P, 0.1, 100.0)
    # vectorize the scalar call
    vrho = np.vectorize(lambda t, p: IAPWS97(T=t, P=p).rho)
    return vrho(T, P_clamped)

def main():
    # Temperature list (°C) and convert to Kelvin
    temps_C = [25, 50, 75, 99]
    temps_K = [t + 273.15 for t in temps_C]

    # Pressure axis: 0.1–100 MPa
    pressures = np.linspace(0.1, 100.0, 500)

    plt.figure(figsize=(10, 6))

    for T_C, T_K in zip(temps_C, temps_K):
        rho_s = rho_simple(T_K, pressures)
        rho_i = rho_iapws(T_K, pressures)
        plt.plot(pressures, rho_s, '--', label=f'Simple {T_C}°C')
        plt.plot(pressures, rho_i, '-',  label=f'IAPWS {T_C}°C')

    plt.title("Water Density vs Pressure at Various Temperatures")
    plt.xlabel("Pressure (MPa)")
    plt.ylabel("Density (kg/m³)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()