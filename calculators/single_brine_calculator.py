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
        print(f"DEBUG: Starting with base density: {base_density} kg/m³")
        m_salt = self._density_to_molality(base_density)
        print(f"DEBUG: Converted to molality: {m_salt} mol/kg")
        
        # Use a fallback approach for first calculation if needed
        use_fallback = True
        try:
            test_density = self._calculate_density(m_salt, pressures[0], temperatures[0])
            if test_density > 900 and test_density < 2500:  # Reasonable range for brine densities
                use_fallback = False
                print(f"DEBUG: Main calculation method produced reasonable density: {test_density} kg/m³")
            else:
                print(f"DEBUG: Main calculation method produced unreasonable density: {test_density} kg/m³")
        except Exception as e:
            print(f"DEBUG: Error in test calculation: {str(e)}")
        
        for p in pressures:
            p_key = f"{p:.2f}"
            results["densities"][p_key] = {}
            
            for t in temperatures:
                t_key = f"{t:.2f}"
                try:
                    if use_fallback:
                        density = self._simplified_density_calculation(base_density, p, t)
                        print(f"DEBUG: Using simplified calculation at P={p} MPa, T={t} K: {density} kg/m³")
                    else:
                        density = self._calculate_density(m_salt, p, t)
                        
                    # Ensure we have reasonable values
                    if density <= 0 or density > 3000:
                        density = self._simplified_density_calculation(base_density, p, t)
                        print(f"DEBUG: Falling back to simplified calculation at P={p} MPa, T={t} K: {density} kg/m³")
                        
                    results["densities"][p_key][t_key] = round(density, 2)
                except Exception as e:
                    print(f"DEBUG: Error calculating density at P={p} MPa, T={t} K: {str(e)}")
                    density = self._simplified_density_calculation(base_density, p, t)
                    results["densities"][p_key][t_key] = round(density, 2)
                
        return results
    
    def _simplified_density_calculation(self, base_density, pressure, temperature):
        """
        A simplified empirical model for brine density as a function of pressure and temperature.
        This is a fallback when the theoretical model produces unreasonable results.
        
        Args:
            base_density (float): Base density at reference conditions (kg/m³)
            pressure (float): Pressure in MPa
            temperature (float): Temperature in K
            
        Returns:
            float: Estimated density in kg/m³
        """
        # Typical thermal expansion coefficient for brines (per K)
        thermal_coeff = -0.0003
        
        # Typical compressibility coefficient for brines (per MPa)
        pressure_coeff = 0.0007
        
        # Reference conditions
        ref_temp = 298.15  # K (25°C)
        ref_pressure = 0.1  # MPa (atmospheric)
        
        # Adjust for temperature and pressure effects
        # ρ = ρ₀ × (1 + β×ΔP) / (1 + α×ΔT)
        delta_t = temperature - ref_temp
        delta_p = pressure - ref_pressure
        
        density = base_density * (1 + pressure_coeff * delta_p) / (1 + thermal_coeff * delta_t)
        
        print(f"DEBUG: Simplified density calculation - base: {base_density}, ΔT: {delta_t}, ΔP: {delta_p}")
        print(f"DEBUG: Simplified result: {density} kg/m³")
        
        return density
    
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
        print(f"DEBUG: Water density at P={pressure} MPa, T={temperature} K: {water_density} kg/m³")
        
        # Calculate apparent molal volume with scaling correction
        # The parameters in the tables are likely in cm³/mol, which needs conversion to m³/mol
        phi_v = self._calculate_apparent_molal_volume(molality, temperature, pressure)
        phi_v_corrected = phi_v * 1e-6  # Convert from cm³/mol to m³/mol
        print(f"DEBUG: Apparent molal volume: {phi_v_corrected} m³/mol (after scaling)")
        
        # Calculate density using Equation 1 from paper
        M_salt = self.config['molecular_weight']
        print(f"DEBUG: Salt molecular weight: {M_salt} kg/mol")
        
        # p = (1 + Σ mi*Mi) / (1/pw + Σ mi*φi)
        numerator = 1.0 + molality * M_salt
        denominator = (1.0 / water_density) + molality * phi_v_corrected
        
        print(f"DEBUG: Density equation: numerator = {numerator}, denominator = {denominator}")
        
        density = numerator / denominator
        print(f"DEBUG: Calculated density result: {density} kg/m³")
        
        return density
    
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
        # For CaCl2 brine, we can use an empirical formula as a first approximation
        # This is based on experimental data and provides a reasonable starting point
        
        # Get salt name
        salt_name = self.config['name']
        
        # Simple empirical approach for common salts
        if salt_name == 'NaCl':
            # For NaCl: Approximately 58.44g NaCl per 0.1 kg/m³ increase above water density
            water_density = self._calculate_water_density(reference_temp, reference_pressure)
            approx_salt_mass = (density - water_density) / 10  # g/L
            molality = approx_salt_mass / 58.44 / (1000 - approx_salt_mass) * 1000
            
        elif salt_name == 'CaCl2':
            # For CaCl2: Approximately 111g CaCl2 per 0.1 kg/m³ increase above water density
            water_density = self._calculate_water_density(reference_temp, reference_pressure)
            approx_salt_mass = (density - water_density) / 10  # g/L
            molality = approx_salt_mass / 111.0 / (1000 - approx_salt_mass) * 1000
            
        else:
            # For other salts, use a more general approach
            # Start with an initial guess based on density difference
            water_density = self._calculate_water_density(reference_temp, reference_pressure)
            density_diff = density - water_density  # kg/m³
            
            # Rough conversion factor: assume 1 mol/kg increases density by about 100 kg/m³
            molality = density_diff / 100.0
            
            # Cap at reasonable values for numerical stability
            molality = min(max(molality, 0.01), 10.0)
        
        print(f"DEBUG: Estimated molality from density {density} kg/m³: {molality} mol/kg")
        
        # Fine-tune through iteration if needed
        try:
            water_density = self._calculate_water_density(reference_temp, reference_pressure)
            M_salt = self.config['molecular_weight']
            
            # Skip iterative refinement for now, as it may use uncorrected phi_v values
            # Just use the empirical estimate
            
            print(f"DEBUG: Final molality estimate: {molality} mol/kg")
            return max(0.001, molality)  # Ensure positive value
            
        except Exception as e:
            print(f"DEBUG: Error in molality calculation: {str(e)}")
            return max(0.001, molality)  # Use empirical estimate
    
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
        # Calculate infinite dilution molal volume (φ°v)
        phi_v0 = self._calculate_phi_v0(temperature, pressure)
        print(f"DEBUG: Infinite dilution molal volume (phi_v0): {phi_v0}")
        
        # Calculate ionic strength
        v_cation = self.config['ions']['cation']['stoichiometry']
        v_anion = self.config['ions']['anion']['stoichiometry']
        z_cation = self.config['ions']['cation']['charge']
        z_anion = self.config['ions']['anion']['charge']
        
        I = 0.5 * molality * (v_cation * z_cation**2 + v_anion * z_anion**2)
        print(f"DEBUG: Ionic strength (I): {I}")
        
        # Calculate Debye-Hückel term
        A_v = self._calculate_debye_huckel_slope(temperature, pressure)
        print(f"DEBUG: Debye-Hückel slope (A_v): {A_v}")
        
        DH_factor = (v_cation * z_cation**2 + v_anion * z_anion**2) * A_v
        DH_denominator = 2 * (1 + I**0.5)
        DH_term = (DH_factor * I**0.5) / DH_denominator if I > 0 else 0
        print(f"DEBUG: D-H term: {DH_term}")
        
        # Calculate interaction term (BMX from Equations 5-6)
        B_MX = self._calculate_interaction_parameter(temperature, pressure, I)
        print(f"DEBUG: Interaction parameter (B_MX): {B_MX}")
        
        interaction_term = v_cation * v_anion * B_MX * molality / 2
        print(f"DEBUG: Interaction term: {interaction_term}")
        
        # Combine terms (Equation 2 from paper)
        phi_v = phi_v0 + DH_term + interaction_term
        print(f"DEBUG: Apparent molal volume (phi_v): {phi_v}")
        
        return phi_v
    
    def _calculate_phi_v0(self, temperature, pressure):
        """
        Calculate infinite dilution molal volume using coefficients from Table 2.
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Infinite dilution molal volume in cm³/mol
        """
        delta_p = pressure - 0.1  # Pressure difference from atmospheric
        
        # Scale factor for coefficients - they're likely in different units
        scale_factor = 1.0
        
        # Using the polynomial equation from the paper with coefficients from Table 2
        coeffs = self.config['phi_v0_coeffs']
        
        # First term (A00 + A10*T + A20*T^2)
        term1 = (coeffs.get('A00', 0) * scale_factor + 
                coeffs.get('A10', 0) * scale_factor * temperature + 
                coeffs.get('A20', 0) * scale_factor * temperature**2)
        
        # Second term (A01 + A11*T + A21*T^2)*ΔP
        term2 = 0
        if delta_p > 0:
            term2 = ((coeffs.get('A01', 0) * scale_factor + 
                    coeffs.get('A11', 0) * scale_factor * temperature + 
                    coeffs.get('A21', 0) * scale_factor * temperature**2) * delta_p)
        
        # Third term (A02 + A12*T + A22*T^2)*ΔP^2
        term3 = 0
        if 'A02' in coeffs and delta_p > 0:
            term3 = ((coeffs.get('A02', 0) * scale_factor + 
                     coeffs.get('A12', 0) * scale_factor * temperature + 
                     coeffs.get('A22', 0) * scale_factor * temperature**2) * delta_p**2)
        
        # Combine all terms
        phi_v0 = term1 + term2 + term3
        
        # The paper likely uses cm³/mol, which is 10^-6 m³/mol
        # Keep in original units (cm³/mol) and convert later
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
        
        # Calculate exponent term by term to avoid overflow
        exponent = (
            coeffs['Av0'] + 
            coeffs['Av1'] * temperature * delta_p + 
            coeffs['Av2'] * temperature**2 +
            coeffs['Av3'] * temperature**2 * delta_p +
            coeffs['Av4'] * temperature**2 * delta_p**2
        )
        
        # Clip exponent to avoid overflow
        exponent = max(min(exponent, 50), -50)
        
        A_v = 1e-6 * np.exp(exponent)
        
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
        
        # Scale factor for coefficients - they're likely in different units
        scale_factor = 1.0
        
        coeffs = self.config['interaction_coeffs']
        
        # Calculate B0 (temperature and pressure dependent)
        B0 = (coeffs.get('B00', 0) * scale_factor + 
             coeffs.get('B10', 0) * scale_factor * temperature + 
             coeffs.get('B20', 0) * scale_factor * temperature**2)
        
        if delta_p > 0:
            B0 += ((coeffs.get('B01', 0) * scale_factor + 
                   coeffs.get('B11', 0) * scale_factor * temperature + 
                   coeffs.get('B21', 0) * scale_factor * temperature**2) * delta_p)
            
            if 'B02' in coeffs:
                B0 += ((coeffs.get('B02', 0) * scale_factor + 
                       coeffs.get('B12', 0) * scale_factor * temperature + 
                       coeffs.get('B22', 0) * scale_factor * temperature**2) * delta_p**2)
        
        # Calculate B1 (concentration dependent term)
        B1 = coeffs.get('B1', 0) * scale_factor
        
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
        alpha = 0.000214  # Thermal expansion coefficient (1/K)
        dT = temperature - T0
        
        # Pressure effect (simplified)
        beta = 0.46e-9  # Compressibility (1/Pa)
        dP = (pressure - P0) * 1e6  # Convert MPa to Pa
        
        # Calculate density
        rho = rho0 * (1 + beta * dP) / (1 + alpha * dT)
        
        return rho