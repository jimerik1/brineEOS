import json
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Load water density data
with open(os.path.join(os.path.dirname(__file__), 'water_density.json')) as f:
    data = json.load(f)

# Parse pressure keys (e.g. "100000.0") as float
available_pressures = sorted([float(p) for p in data['densities'].keys()])
available_temperatures = sorted({
    float(t) for p in data['densities'].values() for t in p.keys()
})

# Convert to target units
pressures = [p / 1e6 for p in available_pressures]  # Pa -> MPa
temperatures = [t + 273.15 for t in available_temperatures]  # °C -> K

# Fill grid with densities
density_grid = np.full((len(pressures), len(temperatures)), np.nan)
for i, p in enumerate(available_pressures):
    p_str = str(p)
    for j, t_C in enumerate(available_temperatures):
        t_str = str(round(t_C, 2))
        density_grid[i][j] = data['densities'][p_str].get(t_str, np.nan)

# Create interpolator
interpolator = RegularGridInterpolator(
    (pressures, temperatures),
    density_grid,
    bounds_error=False,
    fill_value=None
)

def lookup_water_density(temperature_K, pressure_MPa):
    """
    Interpolates water density from the table.

    Args:
        temperature_K (float): Temperature in Kelvin
        pressure_MPa (float): Pressure in MPa

    Returns:
        float: Water density in kg/m³
    """
    return float(interpolator((pressure_MPa, temperature_K)))