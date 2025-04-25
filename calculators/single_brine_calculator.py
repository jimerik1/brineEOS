# calculators/single_brine_calculator.py
import numpy as np
from .brine_utils import (
    calculate_water_density, 
    calculate_debye_huckel_slope, 
    calculate_phi_v0,
    calculate_interaction_parameter
)
from scipy.optimize import brentq

class SingleBrineCalculator:
    def __init__(self, salt_config):
        """
        Initialize the calculator with salt-specific parameters.
        
        Args:
            salt_config (dict): Configuration parameters for the salt
        """
        self.config = salt_config
    
    def calculate(self, base_density, pressure_interval, pressure_resolution, 
                 temperature_interval, temperature_resolution):
        """
        Calculate brine densities for a range of pressures and temperatures.
        
        Args:
            base_density (float): Base density of the brine in kg/m³
            pressure_interval (list): [min_pressure, max_pressure] in MPa
            pressure_resolution (float): Pressure step size in MPa
            temperature_interval (list): [min_temperature, max_temperature] in K
            temperature_resolution (float): Temperature step size in K
            
        Returns:
            dict: Calculated densities at each pressure and temperature point
        """
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
                "brine_type": self.config['name'],
                "base_density": base_density,
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
        
        # Convert base density to molality
        m_salt = self._density_to_molality(base_density)
        print(f"DEBUG: Input base density {base_density} kg/m³ converted to molality {m_salt:.6f} mol/kg")
        
        # Validate that our conversion works by calculating density at reference conditions
        reference_density = self._calculate_density(m_salt, 0.1, 298.15)
        print(f"DEBUG: Calculated density at reference conditions: {reference_density:.2f} kg/m³ (should be close to {base_density:.2f})")
        
        for p in pressures:
            p_key = f"{p:.2f}"
            results["densities"][p_key] = {}
            
            for t in temperatures:
                t_key = f"{t:.2f}"
                density = self._calculate_density(m_salt, p, t)
                results["densities"][p_key][t_key] = round(density, 2)
                
        return results
    
    def _calculate_density(self, molality, pressure, temperature):
        """
        Calculate density at a specific pressure and temperature using Equation 1.
        
        Args:
            molality (float): Salt molality in mol/kg
            pressure (float): Pressure in MPa
            temperature (float): Temperature in K
            
        Returns:
            float: Density in kg/m³
        """
        print(f"\nDEBUG: Calculating density for {self.config['name']} at {molality:.6f} mol/kg, P={pressure}MPa, T={temperature}K")
        
        
        # Calculate water density at given temperature and pressure
        water_density = calculate_water_density(temperature, pressure)
        
        
        # Calculate apparent molal volume
        phi_v = self._calculate_apparent_molal_volume(molality, temperature, pressure)
        print(f"DEBUG: Apparent molal volume (ϕ_v) before unit conversion: {phi_v:.6e} cm³/mol")
        
                
        # Calculate density using Equation 1 from paper
        M_salt = self.config['molecular_weight']
        print(f"DEBUG: Salt molecular weight: {M_salt} kg/mol")
        
        # p = (1 + Σ mi*Mi) / (1/pw + Σ mi*φi)
        numerator = 1.0 + molality * M_salt
        denominator = (1.0 / water_density) + molality * phi_v
        
        print(f"DEBUG: Density equation: numerator = {numerator:.6f}, denominator = {denominator:.6e}")
        
        density = numerator / denominator
        print(f"DEBUG: Calculated density result: {density:.2f} kg/m³")
        
        return density
    
    def _density_to_molality(self, rho_target,
                            T_ref=298.15, P_ref=0.1):

        def f(m):
            return self._calculate_density(m, P_ref, T_ref) - rho_target

        # lower and upper bounds (mol/kg)
        m_lo, m_hi = 1e-4, 10.0          # extend if very heavy brine
        m_root = brentq(f, m_lo, m_hi, xtol=1e-8, rtol=1e-6, maxiter=50)
        return m_root
    
    def _calculate_apparent_molal_volume(self, molality, temperature, pressure):
        """
        Calculate apparent molal volume using the Brønsted-Guggenheim model (Equation 2).
        
        Args:
            molality (float): Salt molality in mol/kg
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Apparent molal volume in cm³/mol (needs to be converted to m³/mol later)
        """
        print(f"DEBUG: Calculating apparent molal volume for {self.config['name']} at {molality:.6f} mol/kg")
        
        # Calculate infinite dilution molal volume (φ°v)
        phi_v0 = calculate_phi_v0(temperature, pressure, self.config['phi_v0_coeffs'])
        
        # Calculate ionic strength
        v_cation = self.config['ions']['cation']['stoichiometry']
        v_anion = self.config['ions']['anion']['stoichiometry']
        z_cation = self.config['ions']['cation']['charge']
        z_anion = self.config['ions']['anion']['charge']
        
        I = 0.5 * molality * (v_cation * z_cation**2 + v_anion * z_anion**2)
        print(f"DEBUG: Ionic strength (I): {I:.4f}")
        
        # Calculate Debye-Hückel term
        A_v = calculate_debye_huckel_slope(temperature, pressure, self.config['debye_huckel_coeffs'])
        
        DH_factor = (v_cation * z_cation**2 + v_anion * z_anion**2) * A_v
        DH_denominator = 2 * (1 + I**0.5)
        DH_term = (DH_factor * I**0.5) / DH_denominator if I > 0 else 0
        print(f"DEBUG: Debye-Hückel term: {DH_term:.6e} cm³/mol")
        
        # Calculate interaction term (BMX from Equations 5-6)
        B_MX = calculate_interaction_parameter(temperature, pressure, I, self.config['interaction_coeffs'])
        
        interaction_term = v_cation * v_anion * B_MX * molality / 2
        print(f"DEBUG: Interaction term: {interaction_term:.6e} cm³/mol")
        
        # Combine terms (Equation 2 from paper)
        phi_v = phi_v0 + DH_term + interaction_term
        print(f"DEBUG: Final apparent molal volume (ϕ_v): {phi_v:.6e} cm³/mol")
        
        return phi_v