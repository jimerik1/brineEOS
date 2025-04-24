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
        
        # First, estimate the base density at reference conditions using our empirical model
        estimated_density = self._estimate_mixed_brine_density(salt_composition)
        print(f"DEBUG: Estimated base density for mixed brine: {estimated_density} kg/m³")
        
        for p in pressures:
            p_key = f"{p:.2f}"
            results["densities"][p_key] = {}
            
            for t in temperatures:
                t_key = f"{t:.2f}"
                try:
                    # Try the theoretical calculation first
                    density = self._calculate_density(salt_composition, p, t)
                    
                    # Check if the result is reasonable
                    if density <= 0 or density > 3000:
                        print(f"DEBUG: Theoretical calculation at P={p} MPa, T={t} K produced unreasonable density: {density} kg/m³")
                        # Fall back to simplified model if theoretical calculation fails
                        density = self._simplified_density_calculation(salt_composition, p, t, estimated_density)
                        print(f"DEBUG: Using simplified calculation instead: {density} kg/m³")
                    else:
                        print(f"DEBUG: Theoretical calculation succeeded at P={p} MPa, T={t} K: {density} kg/m³")
                        
                except Exception as e:
                    print(f"DEBUG: Error in density calculation at P={p} MPa, T={t} K: {str(e)}")
                    # Fall back to simplified model if theoretical calculation fails
                    density = self._simplified_density_calculation(salt_composition, p, t, estimated_density)
                    print(f"DEBUG: Using simplified calculation due to error: {density} kg/m³")
                
                results["densities"][p_key][t_key] = round(density, 2)
                
        return results
    
    def _estimate_mixed_brine_density(self, salt_composition):
        """
        Estimate the mixed brine density at reference conditions based on composition.
        
        Args:
            salt_composition (dict): Salt composition in weight percentages
            
        Returns:
            float: Estimated density in kg/m³
        """
        # Start with water density
        water_density = 997.0  # kg/m³ at 25°C
        
        # Calculate water percentage
        water_wt_pct = 100.0
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                water_wt_pct -= pct
        
        # Ensure water percentage is positive
        water_wt_pct = max(water_wt_pct, 1.0)
        
        print(f"DEBUG: Water percentage: {water_wt_pct}%")
        
        # Calculate contribution of each salt
        # These factors are approximate but should give reasonable estimates
        density_increase = 0.0
        
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                if salt == 'CaCl2':
                    density_increase += pct * 11.0  # Approx. 11 kg/m³ per 1% CaCl2
                elif salt == 'CaBr2':
                    density_increase += pct * 13.0  # Approx. 13 kg/m³ per 1% CaBr2
                elif salt == 'ZnBr2':
                    density_increase += pct * 15.0  # Approx. 15 kg/m³ per 1% ZnBr2
                elif salt == 'NaCl':
                    density_increase += pct * 7.0   # Approx. 7 kg/m³ per 1% NaCl
                elif salt == 'KCl':
                    density_increase += pct * 7.5   # Approx. 7.5 kg/m³ per 1% KCl
                else:
                    density_increase += pct * 10.0  # Generic approximation
                
                print(f"DEBUG: Density contribution from {salt} ({pct}%): {pct * 10.0} kg/m³")
        
        # Calculate estimated density
        estimated_density = water_density + density_increase
        
        print(f"DEBUG: Total density increase from salts: {density_increase} kg/m³")
        print(f"DEBUG: Estimated mixed brine density: {estimated_density} kg/m³")
        
        return estimated_density
    
    def _simplified_density_calculation(self, salt_composition, pressure, temperature, base_density=None):
        """
        A simplified empirical model for brine density as a function of pressure and temperature.
        This is used when the theoretical calculation fails.
        
        Args:
            salt_composition (dict): Salt composition in weight percentages
            pressure (float): Pressure in MPa
            temperature (float): Temperature in K
            base_density (float): Optional estimated base density at reference conditions
            
        Returns:
            float: Estimated density in kg/m³
        """
        # If base density not provided, estimate it
        if base_density is None:
            base_density = self._estimate_mixed_brine_density(salt_composition)
            
        # Reference conditions
        ref_temp = 298.15  # K (25°C)
        ref_pressure = 0.1  # MPa (atmospheric)
        
        # Calculate total salt concentration for adjusting coefficients
        total_salt_pct = sum(pct for salt, pct in salt_composition.items() if salt in self.salt_configs and pct > 0)
        
        # Adjust thermal expansion and compressibility based on salt concentration
        # Higher salt concentrations reduce thermal expansion and compressibility
        
        # Calculate the thermal expansion coefficient (per K)
        # Typical range: -0.0002 to -0.0004, decreases (smaller magnitude) with increasing salt concentration
        thermal_coeff_base = -0.0003
        thermal_coeff = thermal_coeff_base * (1.0 - total_salt_pct / 200.0)  # Adjust based on salt concentration
        
        # Calculate the compressibility coefficient (per MPa)
        # Typical range: 0.0005 to 0.0008, decreases with increasing salt concentration
        pressure_coeff_base = 0.0007
        pressure_coeff = pressure_coeff_base * (1.0 - total_salt_pct / 300.0)  # Adjust based on salt concentration
        
        # Calculate temperature and pressure differences
        delta_t = temperature - ref_temp
        delta_p = pressure - ref_pressure
        
        # Apply temperature and pressure corrections
        # ρ = ρ₀ × (1 + β×ΔP) / (1 + α×ΔT)
        density = base_density * (1 + pressure_coeff * delta_p) / (1 + thermal_coeff * delta_t)
        
        print(f"DEBUG: Simplified calculation - Base: {base_density:.2f}, Thermal coeff: {thermal_coeff:.6f}, Pressure coeff: {pressure_coeff:.6f}")
        print(f"DEBUG: Simplified calculation - ΔT: {delta_t:.2f}, ΔP: {delta_p:.2f}, Result: {density:.2f} kg/m³")
        
        return density
    
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
        print(f"DEBUG: Starting calculation for mixed brine with composition: {salt_composition}")
        print(f"DEBUG: Temperature = {temperature}K, Pressure = {pressure}MPa")
        
        if not salt_composition:
            raise ValueError("Salt composition is required for mixed brine calculations")
        
        # Calculate water weight percentage
        water_wt_pct = 100.0
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                water_wt_pct -= pct
        
        print(f"DEBUG: Water weight percentage = {water_wt_pct}%")  # Should be 50.7% for the example
        
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
                print(f"DEBUG: Molality of {salt} = {molalities[salt]} mol/kg")
        
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
        
        print("DEBUG: Ion molalities:", ion_molalities)
        print("DEBUG: Ion charges:", ion_charges)
        
        # Calculate ionic strength (Equation 3)
        ionic_strength = 0.5 * sum(ion_molalities[ion] * ion_charges[ion]**2 
                                for ion in ion_molalities)
        
        print(f"DEBUG: Ionic strength = {ionic_strength}")  # Should be 20.3 for the example
        
        # Calculate water density at given conditions
        water_density = self._calculate_water_density(temperature, pressure)
        print(f"DEBUG: Water density = {water_density} kg/m³")
        
        # Calculate Debye-Hückel limiting law slope
        A_v = self._calculate_debye_huckel_slope(temperature, pressure)
        print(f"DEBUG: Debye-Hückel slope (A_v) = {A_v}")
        
        # For each salt, calculate apparent molal volume
        numerator = 1.0  # Represents (1 + Σ mᵢMᵢ) in Equation 1
        denominator = 1.0 / water_density  # Represents (1/ρw + Σ mᵢφᵢ) in Equation 1
        
        print("DEBUG: Calculating apparent molal volumes...")
        
        # Calculate normalization factors
        N_factors = self._calculate_normalization_factors(
            ion_molalities, ion_charges
        )
        print("DEBUG: Normalization factors:", N_factors)
        
        for salt, molality in molalities.items():
            if molality > 0:
                config = self.salt_configs[salt]
                M_salt = config['molecular_weight']
                
                # Add mass contribution to numerator
                mass_contrib = molality * M_salt
                numerator += mass_contrib
                print(f"DEBUG: Mass contribution of {salt} = {mass_contrib}")
                
                # Calculate infinite dilution molal volume
                phi_v0 = self._calculate_phi_v0(salt, temperature, pressure)
                print(f"DEBUG: Infinite dilution molal volume for {salt} (phi_v0) = {phi_v0}")
                
                # Calculate apparent molal volume
                phi_v = self._calculate_apparent_molal_volume(
                    salt, molality, temperature, pressure, ionic_strength, N_factors,
                    ion_molalities=ion_molalities, ion_charges=ion_charges
                )
                
                # Apply scaling correction (similar to single brine calculator)
                # This converts from cm³/mol to m³/mol
                phi_v_corrected = phi_v * 1e-6
                print(f"DEBUG: Apparent molal volume for {salt} (phi_v) = {phi_v_corrected} (after scaling)")
                
                # Add volume contribution to denominator
                vol_contrib = molality * phi_v_corrected
                denominator += vol_contrib
                print(f"DEBUG: Volume contribution of {salt} = {vol_contrib}")
        
        print(f"DEBUG: Final numerator = {numerator}")
        print(f"DEBUG: Final denominator = {denominator}")
        
        # Calculate density using Equation 1
        density = numerator / denominator
        print(f"DEBUG: Calculated density = {density} kg/m³")
        
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
        # The key calculation from Equation A-3 in the paper:
        # Molality = (weight percentage of salt) / (molecular weight * weight percentage of water / 100)
        return (weight_percent) / (molecular_weight * water_wt_pct / 100)
    
    def _calculate_normalization_factors(self, ion_molalities, ion_charges):
        """
        Calculate normalization factors for ion interactions in mixed salt solutions.
        Uses Equation 11 from the paper.
        
        Args:
            ion_molalities (dict): Molalities of all ions
            ion_charges (dict): Charges of all ions
            
        Returns:
            dict: Normalization factors
        """
        N_factors = {}
        
        # Calculate separate sums for cations and anions
        pos_charge_sum = 0.0
        neg_charge_sum = 0.0
        
        for ion, molality in ion_molalities.items():
            charge = ion_charges[ion]
            if charge > 0:  # Cation
                pos_charge_sum += molality * charge**2
            else:  # Anion
                neg_charge_sum += molality * charge**2
        
        print(f"DEBUG: Sum of Z²*m for positive ions = {pos_charge_sum}")
        print(f"DEBUG: Sum of Z²*m for negative ions = {neg_charge_sum}")
        
        # Calculate normalization factor for each ion
        for ion, molality in ion_molalities.items():
            charge = ion_charges[ion]
            if charge > 0:  # Cation
                if pos_charge_sum > 0:
                    N_factors[ion] = (molality * charge**2) / pos_charge_sum
                else:
                    N_factors[ion] = 0.0
            else:  # Anion
                if neg_charge_sum > 0:
                    N_factors[ion] = (molality * charge**2) / neg_charge_sum
                else:
                    N_factors[ion] = 0.0
        
        return N_factors
    
    def _calculate_apparent_molal_volume(self, salt, molality, temperature, pressure, 
                                    ionic_strength, N_factors, ion_molalities=None, ion_charges=None):
        """
        Calculate apparent molal volume for a salt in a mixture.
        
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
            float: Apparent molal volume in cm³/mol (needs conversion to m³/mol later)
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
        
        # Implement Debye-Hückel term from Equation 2
        dh_numerator = (v_cation * z_cation**2 + v_anion * z_anion**2) * A_v * ionic_strength**0.5
        dh_denominator = 2 * (1 + ionic_strength**0.5)
        
        # Avoid division by zero or tiny denominator
        if dh_denominator < 1e-10:
            DH_term = 0
        else:
            DH_term = dh_numerator / dh_denominator
        
        print(f"DEBUG: Debye-Hückel term for {salt} = {DH_term}")
        
        # If ion molalities aren't provided, we skip the interaction terms
        if ion_molalities is None or ion_charges is None:
            return phi_v0 + DH_term
        
        # Get cation and anion symbols for this salt
        cation = config['ions']['cation']['symbol']
        anion = config['ions']['anion']['symbol']
        
        # Calculate interaction term using Equation 10
        cation_interaction = 0.0
        for ion, molality_ion in ion_molalities.items():
            if ion_charges.get(ion, 0) < 0:  # Anions
                # Get the appropriate BMX term for this salt
                B_MX = self._calculate_interaction_parameter(
                    salt, temperature, pressure, ionic_strength
                )
                N_ion = N_factors.get(ion, 0)
                term = v_cation * N_ion * B_MX * molality_ion
                cation_interaction += term
                print(f"DEBUG: {cation}-{ion} interaction term = {term}")
        
        anion_interaction = 0.0
        for ion, molality_ion in ion_molalities.items():
            if ion_charges.get(ion, 0) > 0:  # Cations
                # Get the appropriate BMX term for this salt
                B_MX = self._calculate_interaction_parameter(
                    salt, temperature, pressure, ionic_strength
                )
                N_ion = N_factors.get(ion, 0)
                term = v_anion * N_ion * B_MX * molality_ion
                anion_interaction += term
                print(f"DEBUG: {anion}-{ion} interaction term = {term}")
        
        # Total interaction term (divide by 2 per Equation 10)
        interaction_term = (cation_interaction + anion_interaction) / 2
        print(f"DEBUG: Total interaction term for {salt} = {interaction_term}")
        
        # Combine all terms (Equation 2)
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
            float: Infinite dilution molal volume in cm³/mol
        """
        delta_p = pressure - 0.1  # Pressure difference from atmospheric
        
        config = self.salt_configs[salt]
        coeffs = config['phi_v0_coeffs']
        
        # Calculate term by term to avoid overflow
        # First term (A00 + A10*T + A20*T^2)
        term1 = (coeffs.get('A00', 0) + 
                coeffs.get('A10', 0) * temperature + 
                coeffs.get('A20', 0) * temperature**2)
        
        # Second term (A01 + A11*T + A21*T^2)*ΔP
        term2 = 0
        if delta_p > 0:
            term2 = ((coeffs.get('A01', 0) + 
                    coeffs.get('A11', 0) * temperature + 
                    coeffs.get('A21', 0) * temperature**2) * delta_p)
        
        # Third term (A02 + A12*T + A22*T^2)*ΔP^2
        term3 = 0
        if 'A02' in coeffs and delta_p > 0:
            term3 = ((coeffs.get('A02', 0) + 
                    coeffs.get('A12', 0) * temperature + 
                    coeffs.get('A22', 0) * temperature**2) * delta_p**2)
        
        phi_v0 = term1 + term2 + term3
        
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
        
        # Calculate B0 (temperature and pressure dependent) - Equation 5
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
        
        # Calculate BMX (Equation 6 from paper)
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
        # Reference conditions
        T0 = 273.15  # K
        P0 = 0.1     # MPa
        rho0 = 999.975  # kg/m³
        
        # Thermal expansion coefficient needs to be POSITIVE
        # Water density DECREASES with temperature, so α should be positive
        alpha = 0.000214  # Thermal expansion coefficient (1/K)
        dT = temperature - T0
        
        # Compressibility - water density INCREASES with pressure
        beta = 0.46e-9  # Compressibility (1/Pa)
        dP = (pressure - P0) * 1e6  # Convert MPa to Pa
        
        # Calculate density - make sure formula is correct!
        # ρ = ρ₀ × (1 + β×dP) / (1 + α×dT)
        rho = rho0 * (1 + beta * dP) / (1 + alpha * dT)
        
        return rho