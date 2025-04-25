# calculators/water_calculator.py
import numpy as np
from utils.coolprop_utils import calculate_water_density_cp

class WaterCalculator:
    def __init__(self):
        """
        Initialize the WaterCalculator.
        """
        pass # No specific configuration needed for pure water

    def calculate(self, pressure_interval, pressure_resolution,
                 temperature_interval, temperature_resolution):
        """
        Calculate pure water densities for a range of pressures and temperatures
        using CoolProp.

        Args:
            pressure_interval (list): [min_pressure, max_pressure] in MPa
            pressure_resolution (float): Pressure step size in MPa
            temperature_interval (list): [min_temperature, max_temperature] in K
            temperature_resolution (float): Temperature step size in K

        Returns:
            dict: Calculated densities at each pressure and temperature point
        """
        if __debug__:
            print(f"DEBUG: Starting pure water calculation.")

        # Generate pressure and temperature points
        pressures = np.arange(
            pressure_interval[0],
            pressure_interval[1] + pressure_resolution,
            pressure_resolution
        )

        temperatures = np.arange(
            temperature_interval[0],
            temperature_interval[1] + temperature_resolution,
            temperature_resolution
        )

        results = {
            "metadata": {
                "fluid": "Water",
                "source": "CoolProp",
                "pressure_points": pressures.tolist(),
                "temperature_points": temperatures.tolist(),
                "units": {
                    "density": "kg/m³",
                    "pressure": "MPa",
                    "temperature": "K"
                }
            },
            "densities": {}
        }

        for p in pressures:
            p_key = f"{p:.2f}"
            results["densities"][p_key] = {}

            for t in temperatures:
                t_key = f"{t:.2f}"
                try:
                    density = calculate_water_density_cp(t, p)
                    results["densities"][p_key][t_key] = round(density, 2)
                except Exception as e:
                    # Handle potential CoolProp errors for specific P, T points
                    print(f"WARNING: CoolProp failed for Water at P={p} MPa, T={t} K: {e}")
                    results["densities"][p_key][t_key] = None # Indicate calculation failed

        return results
