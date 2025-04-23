# calculators/single_brine_calculator.py
import numpy as np

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
        
        # Convert base density to salt concentration (molality)
        m_salt = self._density_to_molality(base_density)
        
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
        # Calculate water density at given temperature and pressure
        water_density = self._calculate_water_density(temperature, pressure)
        
        # Calculate apparent molal volume
        phi_v = self._calculate_apparent_molal_volume(molality, temperature, pressure)
        
        # Calculate density using Equation 1 from paper
        M_salt = self.config['molecular_weight']
        
        # p = (1 + Σ mi*Mi) / (1/pw + Σ mi*φi)
        numerator = 1.0 + molality * M_salt
        denominator = (1.0 / water_density) + molality * phi_v
        
        return numerator / denominator
    
    def _density_to_molality(self, density, reference_temp=298.15, reference_pressure=0.1):
        """
        Convert density to molality by iterative solution of Equation 1.
        
        Args:
            density (float): Density in kg/m³
            reference_temp (float): Reference temperature in K
            reference_pressure (float): Reference pressure in MPa
            
        Returns:
            float: Molality in mol/kg
        """
        # Calculate water density at reference conditions
        water_density = self._calculate_water_density(reference_temp, reference_pressure)
        
        # Iterative approach to find molality
        m = 0.1  # Initial guess
        M_salt = self.config['molecular_weight']
        
        for _ in range(10):  # Usually converges in a few iterations
            phi_v = self._calculate_apparent_molal_volume(m, reference_temp, reference_pressure)
            
            # From Equation 1: p = (1 + m*M) / (1/pw + m*φ)
            # Rearranged: m = (p/pw - 1) / (p*φ - M)
            m_new = (density/water_density - 1) / (density*phi_v/1000 - M_salt)
            
            if abs(m - m_new) < 0.001:
                m = m_new
                break
            
            # Damped update for stability
            m = 0.7 * m + 0.3 * m_new
        
        return max(0.0, m)
    
    def _calculate_apparent_molal_volume(self, molality, temperature, pressure):
        """
        Calculate apparent molal volume using the Brønsted-Guggenheim model (Equation 2).
        
        Args:
            molality (float): Salt molality in mol/kg
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Apparent molal volume in m³/mol
        """
        # Calculate infinite dilution molal volume (φ°v)
        phi_v0 = self._calculate_phi_v0(temperature, pressure)
        
        # Calculate ionic strength
        v_cation = self.config['ions']['cation']['stoichiometry']
        v_anion = self.config['ions']['anion']['stoichiometry']
        z_cation = self.config['ions']['cation']['charge']
        z_anion = self.config['ions']['anion']['charge']
        
        I = 0.5 * molality * (v_cation * z_cation**2 + v_anion * z_anion**2)
        
        # Calculate Debye-Hückel term
        A_v = self._calculate_debye_huckel_slope(temperature, pressure)
        
        numerator = (v_cation * z_cation**2 + v_anion * z_anion**2) * A_v * I**0.5
        denominator = 2 * (1 + I**0.5)
        DH_term = numerator / denominator
        
        # Calculate interaction term (BMX from Equations 5-6)
        B_MX = self._calculate_interaction_parameter(temperature, pressure, I)
        interaction_term = v_cation * v_anion * B_MX * molality / 2
        
        # Combine terms (Equation 2 from paper)
        phi_v = phi_v0 + DH_term + interaction_term
        
        return phi_v
    
    def _calculate_phi_v0(self, temperature, pressure):
        """
        Calculate infinite dilution molal volume using coefficients from Table 2.
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Infinite dilution molal volume in m³/mol
        """
        delta_p = pressure - 0.1  # Pressure difference from atmospheric
        
        # Using the polynomial equation from the paper with coefficients from Table 2
        coeffs = self.config['phi_v0_coeffs']
        
        phi_v0 = (coeffs['A00'] + 
                 coeffs.get('A10', 0) * temperature + 
                 coeffs.get('A20', 0) * temperature**2)
        
        if delta_p > 0:
            phi_v0 += ((coeffs.get('A01', 0) + 
                       coeffs.get('A11', 0) * temperature + 
                       coeffs.get('A21', 0) * temperature**2) * delta_p)
            
            if 'A02' in coeffs and delta_p > 0:
                phi_v0 += ((coeffs.get('A02', 0) + 
                           coeffs.get('A12', 0) * temperature + 
                           coeffs.get('A22', 0) * temperature**2) * delta_p**2)
        
        return phi_v0
    
    def _calculate_debye_huckel_slope(self, temperature, pressure):
        """
        Calculate Debye-Hückel limiting law slope using Equation 4.
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Debye-Hückel limiting law slope
        """
        delta_p = pressure - 0.1
        
        coeffs = self.config['debye_huckel_coeffs']
        
        A_v = 1e-6 * np.exp(
            coeffs['Av0'] + 
            coeffs['Av1'] * temperature * delta_p + 
            coeffs['Av2'] * temperature**2 +
            coeffs['Av3'] * temperature**2 * delta_p +
            coeffs['Av4'] * temperature**2 * delta_p**2
        )
        
        return A_v
    
    def _calculate_interaction_parameter(self, temperature, pressure, ionic_strength):
        """
        Calculate interaction parameter BMX using coefficients from Table 3.
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            ionic_strength (float): Ionic strength
            
        Returns:
            float: Interaction parameter
        """
        delta_p = pressure - 0.1
        
        coeffs = self.config['interaction_coeffs']
        
        # Calculate B0 (temperature and pressure dependent)
        B0 = (coeffs.get('B00', 0) + 
             coeffs.get('B10', 0) * temperature + 
             coeffs.get('B20', 0) * temperature**2)
        
        if delta_p > 0:
            B0 += ((coeffs.get('B01', 0) + 
                   coeffs.get('B11', 0) * temperature + 
                   coeffs.get('B21', 0) * temperature**2) * delta_p)
            
            if 'B02' in coeffs and delta_p > 0:
                B0 += ((coeffs.get('B02', 0) + 
                       coeffs.get('B12', 0) * temperature + 
                       coeffs.get('B22', 0) * temperature**2) * delta_p**2)
        
        # Calculate B1 (concentration dependent term)
        B1 = coeffs.get('B1', 0)
        
        # Calculate BMX (Equation 5-6 from paper)
        B_MX = B0 + B1 * ionic_strength
        
        return B_MX
    
    def _calculate_water_density(self, temperature, pressure):
        """
        Calculate the density of pure water at given T and P.
        Uses a simplified model valid for 273.15K ≤ T ≤ 647K and P ≤ 100MPa
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Water density in kg/m³
        """
        # Reference conditions
        T0 = 273.15  # K
        P0 = 0.1     # MPa
        rho0 = 999.975  # kg/m³
        
        # Temperature effect (simplified)
        alpha = -0.000214  # Thermal expansion coefficient (1/K)
        dT = temperature - T0
        
        # Pressure effect (simplified)
        beta = 0.46e-9  # Compressibility (1/Pa)
        dP = (pressure - P0) * 1e6  # Convert MPa to Pa
        
        # Calculate density
        rho = rho0 * (1 + beta * dP) / (1 + alpha * dT)
        
        return rho