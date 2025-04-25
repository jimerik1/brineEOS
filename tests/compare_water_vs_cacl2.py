#!/usr/bin/env python3
"""
Compare pure-water density and a CaCl2 brine (ρ0 = 1300 kg/m3)
from 0-100 MPa at three temperatures: 298.15 K, 373.15 K, 423.15 K
(25 °C, 100 °C, 150 °C).

Outputs a PNG named 'water_vs_cacl2.png'.
"""

import requests
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

# ------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------
BASE_URL   = "http://localhost:5099/api/v1"
T_KELVIN   = [298.15, 373.15, 423.15]     # 25, 100, 150 °C
P_MIN, P_MAX = 0.5, 100.0                # MPa
P_STEP       = 1.0

HEADERS = {"Content-Type": "application/json"}

# ------------------------------------------------------------------
# helper to hit endpoint and unpack densities
# ------------------------------------------------------------------
def query(endpoint: str, payload: dict) -> dict:
    url = f"{BASE_URL}{endpoint}"
    r = requests.post(url, json=payload, headers=HEADERS, timeout=30)
    try:
        r.raise_for_status()
    except requests.HTTPError as exc:
        print(f"[{endpoint}] {exc}\n{r.text}")
        sys.exit(1)
    return r.json()

def extract_line(resp: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return pressures and the single-temperature density values."""
    pressures = np.array([float(p) for p in resp["densities"].keys()])
    # only one temperature in each payload → grab first key inside dict
    t_key = next(iter(resp["densities"][next(iter(resp['densities']))]))
    densities = np.array([resp["densities"][f"{p:.2f}"][t_key] for p in pressures])
    return pressures, densities

# ------------------------------------------------------------------
# gather data
# ------------------------------------------------------------------
pressures = np.arange(P_MIN, P_MAX + P_STEP, P_STEP)

water_payload = {
    "pressure_interval": [P_MIN, P_MAX],
    "pressure_resolution": P_STEP,
    # temperature loop will replace these two values each time
    "temperature_interval": [0, 0],
    "temperature_resolution": 1.0e-4
}

brine_payload_template = {
    "brine_type": "CaCl2",
    "base_density": 1300.0,
    "pressure_interval": [P_MIN, P_MAX],
    "pressure_resolution": P_STEP,
    "temperature_interval": [0, 0],
    "temperature_resolution": 1.0e-4
}

series = []   # list of (label, pressure_array, density_array, style)

for T in T_KELVIN:
    # --- water ---
    water_payload["temperature_interval"] = [T, T + 1.0e-4]
    resp_w = query("/calculate_water_density", water_payload)
    p, rho_w = extract_line(resp_w)
    series.append((f"Water {T-273.15:.0f}°C", p, rho_w, {"color": "tab:blue"}))

    # --- CaCl2 brine ---
    brine_payload = brine_payload_template.copy()
    brine_payload["temperature_interval"] = [T, T + 1.0e-4]
    resp_b = query("/calculate_density", brine_payload)
    p, rho_b = extract_line(resp_b)
    series.append((f"CaCl₂ 1300 kg/m³ {T-273.15:.0f}°C", p, rho_b,
                   {"color": "tab:orange", "linestyle": "--"}))

# ------------------------------------------------------------------
# plotting
# ------------------------------------------------------------------
plt.figure(figsize=(8, 6))

for label, p, rho, style in series:
    plt.plot(p, rho, label=label, **style)

plt.title("Water vs CaCl₂ Brine Density  (0–100 MPa)")
plt.xlabel("Pressure (MPa)")
plt.ylabel("Density (kg m$^{-3}$)")
plt.grid(alpha=0.3)
plt.legend(fontsize=9)
plt.tight_layout()

out = Path("water_vs_cacl2.png")
plt.savefig(out, dpi=300)
print(f"Plot saved → {out.resolve()}")

plt.show()