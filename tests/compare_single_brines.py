#!/usr/bin/env python3
"""
Plot densities of every single-salt brine (NaCl, KCl, CaCl2, CaBr2, ZnBr2, ZnCl2)
from 0–100 MPa at 25 °C, 100 °C, 150 °C, using points only (no connecting lines).
"""

import requests
import matplotlib.pyplot as plt
import numpy as np
import itertools
from pathlib import Path

# --------------------------- configuration ---------------------------------
BASE_URL        = "http://localhost:5099/api/v1/calculate_density"
SALTS           = ["NaCl", "KCl", "CaCl2", "CaBr2", "ZnBr2", "ZnCl2"]
T_K             = [298.15, 373.15, 423.15]          # 25 °C, 100 °C, 150 °C
P_MIN, P_MAX    = 0.1, 100.0
P_STEP          = 1.0
BASE_DENSITY    = 1200.0                            # kg m⁻³ for all salts
CMAP            = plt.get_cmap("tab20")
HEADERS         = {"Content-Type": "application/json"}

# --------------------------- helpers ---------------------------------------
def query(payload: dict) -> dict:
    r = requests.post(BASE_URL, json=payload, headers=HEADERS, timeout=30)
    r.raise_for_status()
    return r.json()

def extract_line(resp: dict) -> tuple[np.ndarray, np.ndarray]:
    pressures = np.array([float(p) for p in resp["densities"].keys()])
    t_key     = next(iter(resp["densities"][next(iter(resp["densities"]))]))
    densities = np.array([resp["densities"][f"{p:.2f}"][t_key] for p in pressures])
    return pressures, densities

# --------------------------- data gathering --------------------------------
series = []   # list of (label, p, rho, style)

color_cycle = itertools.cycle(CMAP.colors)

for T in T_K:
    for salt in SALTS:
        payload = {
            "brine_type": salt,
            "base_density": BASE_DENSITY,
            "pressure_interval": [P_MIN, P_MAX],
            "pressure_resolution": P_STEP,
            "temperature_interval": [T, T + 1e-4],
            "temperature_resolution": 1e-4
        }
        resp = query(payload)
        p, rho = extract_line(resp)

        colour = next(color_cycle)
        marker = "o"           # use filled circles for every point
        size   = 25            # marker size

        label = f"{salt}  {T-273.15:.0f} °C"
        series.append((label, p, rho,
                       dict(color=colour, marker=marker, s=size)))

# --------------------------- plotting --------------------------------------
plt.figure(figsize=(10, 7))

for label, p, rho, style in series:
    # scatter: points only, no lines
    plt.scatter(p, rho, label=label, **style, alpha=0.8)

plt.title("Single-Salt Brine Densities (base ρ₀ = 1200 kg m⁻³)\n"
          "0 – 100 MPa at 25 °C, 100 °C, 150 °C")
plt.xlabel("Pressure (MPa)")
plt.ylabel("Density (kg m$^{-3}$)")
plt.grid(alpha=0.3)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()

out = Path("all_single_brines_points.png")
plt.savefig(out, dpi=300)
print(f"Plot saved → {out.resolve()}")
plt.show()