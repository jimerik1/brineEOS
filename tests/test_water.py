#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import requests

API_URL = "http://localhost:5099/api/v1/calculate_water_density"

# Temperatures to plot (in Kelvin)
temps_K = [295.0, 366.0, 447.0]

# Pressure grid: 0.1–150.0 MPa in 0.1 MPa steps (→ 1–1500 bar)
min_p, max_p, dp = 1.0, 150.0, 0.1

def fetch_density_at_T(T_K):
    # Build a small 0.1 K interval *above* T_K when possible,
    # but if T_K is at the 447 K limit, step *below* instead.
    if T_K < 447.0:
        t_low, t_high = T_K, T_K + 0.1
    else:
        t_low, t_high = T_K - 0.1, T_K

    payload = {
        "pressure_interval": [min_p, max_p],
        "pressure_resolution": dp,
        "temperature_interval": [t_low, t_high],
        "temperature_resolution": round(abs(t_high - t_low), 3)
    }

    resp = requests.post(API_URL, json=payload)
    # if you still get 400, uncomment next line to see the details:
    # print(resp.status_code, resp.text)
    resp.raise_for_status()

    data = resp.json()
    pressures = data['metadata']['pressure_points']  # in MPa
    t_key = f"{T_K:.2f}"

    densities = [
        data['densities'][f"{p:.2f}"].get(t_key)
        for p in pressures
    ]
    return pressures, densities

def main():
    plt.figure(figsize=(10, 6))

    for T_K in temps_K:
        ps, ds = fetch_density_at_T(T_K)
        ps_bar = np.array(ps) * 10.0  # MPa → bar
        plt.plot(ps_bar, ds, label=f"{int(T_K)} K", linewidth=1.5)

    plt.title("Water Density vs Pressure at 295 K, 366 K, 447 K")
    plt.xlabel("Pressure (bar)")
    plt.ylabel("Density (kg/m³)")
    plt.legend(title="Temperature")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()