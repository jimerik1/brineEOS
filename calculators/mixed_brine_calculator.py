# calculators/mixed_brine_calculator.py
import numpy as np
from .brine_utils import (
    calculate_water_density, 
    calculate_debye_huckel_slope, 
    calculate_phi_v0,
    calculate_interaction_parameter
)

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
        print(f"DEBUG: Starting mixed brine calculation with composition: {salt_composition}")
        
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
        print(f"\nDEBUG: Calculating density for mixed brine at P={pressure}MPa, T={temperature}K")
        print(f"DEBUG: Salt composition: {salt_composition}")
        
        if not salt_composition:
            raise ValueError("Salt composition is required for mixed brine calculations")
        
        # Calculate water weight percentage
        water_wt_pct = 100.0
        for salt, pct in salt_composition.items():
            if salt in self.salt_configs and pct > 0:
                water_wt_pct -= pct
        
        print(f"DEBUG: Water weight percentage: {water_wt_pct:.2f}%")
        
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
                print(f"DEBUG: Molality of {salt}: {molalities[salt]:.6f} mol/kg")
        
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
        
        print(f"DEBUG: Ion molalities: {ion_molalities}")
        print(f"DEBUG: Ion charges: {ion_charges}")
        
        # Calculate ionic strength (Equation 3)
        ionic_strength = 0.5 * sum(ion_molalities[ion] * ion_charges[ion]**2 
                                for ion in ion_molalities)
        print(f"DEBUG: Ionic strength: {ionic_strength:.4f}")
        
        # Calculate water density at given conditions
        water_density = calculate_water_density(temperature, pressure)
        
        # Calculate normalization factors
        N_factors = self._calculate_normalization_factors(
            ion_molalities, ion_charges
        )
        print(f"DEBUG: Normalization factors: {N_factors}")
        
        # For each salt, calculate apparent molal volume
        numerator = 1.0  # Represents (1 + Σ mᵢMᵢ) in Equation 1
        denominator = 1.0 / water_density  # Represents (1/ρw + Σ mᵢφᵢ) in Equation 1
        
        for salt, molality in molalities.items():
            if molality > 0:
                config = self.salt_configs[salt]
                M_salt = config['molecular_weight']
                
                # Add mass contribution to numerator
                mass_contrib = molality * M_salt
                numerator += mass_contrib
                print(f"DEBUG: Mass contribution of {salt}: {mass_contrib:.6f}")
                
                # Calculate apparent molal volume
                phi_v = self._calculate_apparent_molal_volume(
                    salt, molality, temperature, pressure, ionic_strength, N_factors,
                    ion_molalities=ion_molalities, ion_charges=ion_charges
                )
                print(f"DEBUG: Apparent molal volume for {salt} before unit conversion: {phi_v:.6e} cm³/mol")
                
                # Convert from cm³/mol to m³/mol
                phi_v_corrected = phi_v * 1e-6
                print(f"DEBUG: Apparent molal volume for {salt} after unit conversion: {phi_v_corrected:.6e} m³/mol")
                
                # Add volume contribution to denominator
                vol_contrib = molality * phi_v_corrected
                denominator += vol_contrib
                print(f"DEBUG: Volume contribution of {salt}: {vol_contrib:.6e}")
        
        print(f"DEBUG: Final numerator: {numerator:.6f}")
        print(f"DEBUG: Final denominator: {denominator:.6e}")
        
        # Calculate density using Equation 1
        density = numerator / denominator
        print(f"DEBUG: Final density: {density:.2f} kg/m³")
        
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
        # The calculation from Equation A-3 in the paper:
        # Molality = (weight percentage of salt) / (molecular weight * weight percentage of water / 100)
        molality = (weight_percent/100) / (molecular_weight * water_wt_pct/100)
        print(f"DEBUG: Converted {weight_percent}% {molecular_weight}kg/mol salt to {molality:.6f} mol/kg")
        return molality
    
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
        
        print(f"DEBUG: Sum of Z²*m for positive ions: {pos_charge_sum:.6f}")
        print(f"DEBUG: Sum of Z²*m for negative ions: {neg_charge_sum:.6f}")
        
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
        print(f"DEBUG: Calculating apparent molal volume for {salt} at {molality:.6f} mol/kg")
        
        config = self.salt_configs[salt]
        
        # Calculate infinite dilution molal volume (φ°v)
        phi_v0 = calculate_phi_v0(temperature, pressure, config['phi_v0_coeffs'])
        
        # Calculate Debye-Hückel term
        A_v = calculate_debye_huckel_slope(temperature, pressure, config['debye_huckel_coeffs'])
        
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
        
        print(f"DEBUG: Debye-Hückel term: {DH_term:.6e} cm³/mol")
        
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
                B_MX = calculate_interaction_parameter(
                    temperature, pressure, ionic_strength, config['interaction_coeffs']
                )
                N_ion = N_factors.get(ion, 0)
                term = v_cation * N_ion * B_MX * molality_ion
                cation_interaction += term
                print(f"DEBUG: {cation}-{ion} interaction term: {term:.6e} cm³/mol")
        
        anion_interaction = 0.0
        for ion, molality_ion in ion_molalities.items():
            if ion_charges.get(ion, 0) > 0:  # Cations
                # Get the appropriate BMX term for this salt
                B_MX = calculate_interaction_parameter(
                    temperature, pressure, ionic_strength, config['interaction_coeffs']
                )
                N_ion = N_factors.get(ion, 0)
                term = v_anion * N_ion * B_MX * molality_ion
                anion_interaction += term
                print(f"DEBUG: {anion}-{ion} interaction term: {term:.6e} cm³/mol")
        
        # Total interaction term (divide by 2 per Equation 10)
        interaction_term = (cation_interaction + anion_interaction) / 2
        print(f"DEBUG: Total interaction term: {interaction_term:.6e} cm³/mol")
        
        # Combine all terms (Equation 2)
        phi_v = phi_v0 + DH_term + interaction_term
        print(f"DEBUG: Final apparent molal volume (ϕ_v): {phi_v:.6e} cm³/mol")
        
        return phi_v