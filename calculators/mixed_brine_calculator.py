# calculators/mixed_brine_calculator.py
import numpy as np

class MixedBrineCalculator:
    def __init__(self, salt_configs):
        """
        Initialize the calculator with configurations for all salt types.
        
        Args:
            salt_configs (dict): Dictionary of salt configurations keyed by salt name
        """
        self.salt_configs = salt_configs
    
    def calculate(self, salt_composition, pressure_interval, pressure_resolution, 
                 temperature_interval, temperature_resolution, base_density=None):
        """
        Calculate mixed brine densities for a range of pressures and temperatures.
        
        Args:
            salt_composition (dict): Composition of the brine in weight percentages
            pressure_interval (list): [min_pressure, max_pressure] in MPa
            pressure_resolution (float): Pressure step size in MPa
            temperature_interval (list): [min_temperature, max_temperature] in K
            temperature_resolution (float): Temperature step size in K
            base_density (float): Optional base density (not used for mixed brines)
            
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
                "brine_type": "mixed",
                "salt_composition": salt_composition,
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
                density = self._calculate_density(salt_composition, p, t)
                results["densities"][p_key][t_key] = round(density, 2)
                
        return results
    
    def _calculate_density(self, salt_composition, pressure, temperature):
        """
        Calculate mixed salt brine density according to Equation 1 and related formulations.
        
        Args:
            salt_composition (dict): Salt composition in weight percentages
            pressure (float): Pressure in MPa
            temperature (float): Temperature in K
            
        Returns:
            float: Density in kg/m³
        """
        if not salt_composition:
            raise ValueError("Salt composition is required for mixed brine calculations")
        
        # Calculate water weight percentage
        water_wt_pct = 100.0
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                water_wt_pct -= pct
        
        if water_wt_pct <= 0:
            raise ValueError("Invalid salt composition: total exceeds 100%")
        
        # Convert weight percentages to molalities
        molalities = {}
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                molalities[salt] = self._wt_pct_to_molality(
                    pct, 
                    self.salt_configs[salt]['molecular_weight'], 
                    water_wt_pct
                )
        
        # Calculate ion molalities and charges
        ion_molalities = {}
        ion_charges = {}
        
        for salt, molality in molalities.items():
            config = self.salt_configs[salt]
            
            # Add cation
            cation = config['ions']['cation']['symbol']
            cation_stoich = config['ions']['cation']['stoichiometry']
            if cation not in ion_molalities:
                ion_molalities[cation] = 0
                ion_charges[cation] = config['ions']['cation']['charge']
            ion_molalities[cation] += molality * cation_stoich
            
            # Add anion
            anion = config['ions']['anion']['symbol']
            anion_stoich = config['ions']['anion']['stoichiometry']
            if anion not in ion_molalities:
                ion_molalities[anion] = 0
                ion_charges[anion] = config['ions']['anion']['charge']
            ion_molalities[anion] += molality * anion_stoich
        
        # Calculate ionic strength (Equation 3)
        ionic_strength = 0.5 * sum(ion_molalities[ion] * ion_charges[ion]**2 
                                for ion in ion_molalities)
        
        # Calculate water density
        water_density = self._calculate_water_density(temperature, pressure)
        
        # For each salt, calculate apparent molal volume
        total_volume_contribution = 0.0
        total_mass_contribution = 0.0
        
        for salt, molality in molalities.items():
            if molality > 0:
                config = self.salt_configs[salt]
                
                # Calculate normalization factors for mixed salt interactions
                N_factors = self._calculate_normalization_factors(
                    salt, ion_molalities, ion_charges
                )
                
                # Calculate apparent molal volume for this salt in the mixture
                # Pass ion_molalities and ion_charges to the method
                phi_v = self._calculate_apparent_molal_volume(
                    salt, molality, temperature, pressure, ionic_strength, N_factors,
                    ion_molalities=ion_molalities, ion_charges=ion_charges
                )
                
                # Add contribution to total volume and mass
                M_salt = config['molecular_weight']
                total_volume_contribution += molality * phi_v
                total_mass_contribution += molality * M_salt
        
        # Calculate density using Equation 1
        density = (1.0 + total_mass_contribution) / ((1.0 / water_density) + total_volume_contribution)
        
        return density
    
    def _wt_pct_to_molality(self, weight_percent, molecular_weight, water_wt_pct):
        """
        Convert weight percentage to molality.
        
        Args:
            weight_percent (float): Weight percentage of salt
            molecular_weight (float): Molecular weight of salt in kg/mol
            water_wt_pct (float): Weight percentage of water
            
        Returns:
            float: Molality in mol/kg
        """
        return (weight_percent * 10.0) / (molecular_weight * water_wt_pct)
    
    def _calculate_normalization_factors(self, salt, ion_molalities, ion_charges):
        """
        Calculate normalization factors for ion interactions in mixed salt solutions.
        Uses Equation 11 from the paper.
        
        Args:
            salt (str): Salt identifier
            ion_molalities (dict): Molalities of all ions
            ion_charges (dict): Charges of all ions
            
        Returns:
            dict: Normalization factors
        """
        N_factors = {}
        
        config = self.salt_configs[salt]
        cation = config['ions']['cation']['symbol']
        anion = config['ions']['anion']['symbol']
        
        # Calculate normalization factors using Equation 11
        # For cations
        total_positive_charge = sum(
            ion_molalities[ion] * ion_charges[ion]**2 
            for ion in ion_molalities if ion_charges[ion] > 0
        )
        
        if total_positive_charge > 0:
            N_factors[cation] = (ion_molalities[cation] * ion_charges[cation]**2) / total_positive_charge
        else:
            N_factors[cation] = 0.0
        
        # For anions
        total_negative_charge = sum(
            ion_molalities[ion] * ion_charges[ion]**2 
            for ion in ion_molalities if ion_charges[ion] < 0
        )
        
        if total_negative_charge > 0:
            N_factors[anion] = (ion_molalities[anion] * ion_charges[anion]**2) / total_negative_charge
        else:
            N_factors[anion] = 0.0
        
        return N_factors
    
    def _calculate_apparent_molal_volume(self, salt, molality, temperature, pressure, 
                                    ionic_strength, N_factors, ion_molalities=None, ion_charges=None):
        """
        Calculate apparent molal volume for a salt in a mixture.
        Uses the same approach as for single salt but with normalization factors.
        
        Args:
            salt (str): Salt identifier
            molality (float): Salt molality in mol/kg
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            ionic_strength (float): Total ionic strength of the solution
            N_factors (dict): Normalization factors for ions
            ion_molalities (dict, optional): Molalities of all ions
            ion_charges (dict, optional): Charges of all ions
            
        Returns:
            float: Apparent molal volume in m³/mol
        """
        config = self.salt_configs[salt]
        
        # Calculate infinite dilution molal volume (φ°v)
        phi_v0 = self._calculate_phi_v0(salt, temperature, pressure)
        
        # Calculate Debye-Hückel term
        A_v = self._calculate_debye_huckel_slope(temperature, pressure)
        
        v_cation = config['ions']['cation']['stoichiometry']
        v_anion = config['ions']['anion']['stoichiometry']
        z_cation = config['ions']['cation']['charge']
        z_anion = config['ions']['anion']['charge']
        
        numerator = (v_cation * z_cation**2 + v_anion * z_anion**2) * A_v * ionic_strength**0.5
        denominator = 2 * (1 + ionic_strength**0.5)
        DH_term = numerator / denominator
        
        # If ion molalities aren't provided, we skip the interaction terms
        if ion_molalities is None or ion_charges is None:
            return phi_v0 + DH_term
        
        # Calculate interaction terms for mixed salts (Equation 10)
        cation = config['ions']['cation']['symbol']
        anion = config['ions']['anion']['symbol']
        
        # Interaction with all anions
        cation_interaction = 0
        for ion, m_ion in ion_molalities.items():
            if ion_charges.get(ion, 0) < 0:  # Only consider anions
                B_MX = self._calculate_interaction_parameter(
                    salt, temperature, pressure, ionic_strength
                )
                N_ion = N_factors.get(ion, 0)
                cation_interaction += v_cation * N_ion * B_MX * m_ion
        
        # Interaction with all cations
        anion_interaction = 0
        for ion, m_ion in ion_molalities.items():
            if ion_charges.get(ion, 0) > 0:  # Only consider cations
                B_MX = self._calculate_interaction_parameter(
                    salt, temperature, pressure, ionic_strength
                )
                N_ion = N_factors.get(ion, 0)
                anion_interaction += v_anion * N_ion * B_MX * m_ion
        
        # Total interaction term
        interaction_term = (cation_interaction + anion_interaction) / 2
        
        # Combine terms
        phi_v = phi_v0 + DH_term + interaction_term
        
        return phi_v
    
    def _calculate_phi_v0(self, salt, temperature, pressure):
        """
        Calculate infinite dilution molal volume for a salt.
        
        Args:
            salt (str): Salt identifier
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Infinite dilution molal volume in m³/mol
        """
        delta_p = pressure - 0.1  # Pressure difference from atmospheric
        
        config = self.salt_configs[salt]
        coeffs = config['phi_v0_coeffs']
        
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
        
        # Use the coefficients from any salt config (they should be the same)
        config = next(iter(self.salt_configs.values()))
        coeffs = config['debye_huckel_coeffs']
        
        A_v = 1e-6 * np.exp(
            coeffs['Av0'] + 
            coeffs['Av1'] * temperature * delta_p + 
            coeffs['Av2'] * temperature**2 +
            coeffs['Av3'] * temperature**2 * delta_p +
            coeffs['Av4'] * temperature**2 * delta_p**2
        )
        
        return A_v
    
    def _calculate_interaction_parameter(self, salt, temperature, pressure, ionic_strength):
        """
        Calculate interaction parameter BMX.
        
        Args:
            salt (str): Salt identifier
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            ionic_strength (float): Ionic strength
            
        Returns:
            float: Interaction parameter
        """
        delta_p = pressure - 0.1
        
        config = self.salt_configs[salt]
        coeffs = config['interaction_coeffs']
        
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
        
        Args:
            temperature (float): Temperature in K
            pressure (float): Pressure in MPa
            
        Returns:
            float: Water density in kg/m³
        """
        # Same implementation as in SingleBrineCalculator
        T0 = 273.15  # K
        P0 = 0.1     # MPa
        rho0 = 999.975  # kg/m³
        
        alpha = -0.0002  # Thermal expansion coefficient (1/K)
        dT = temperature - T0
        
        beta = 0.46e-9  # Compressibility (1/Pa)
        dP = (pressure - P0) * 1e6  # Convert MPa to Pa
        
        rho = rho0 * (1 + beta * dP) / (1 + alpha * dT)
        
        return rho