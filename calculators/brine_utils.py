# calculators/brine_utils.py
import numpy as np
import math
from iapws import IAPWS97
def calculate_water_density(temperature, pressure):
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
    
    # Temperature effect
    alpha = 0.000214  # Thermal expansion coefficient (1/K)
    dT = temperature - T0
    
    # Pressure effect
    beta = 0.46e-9  # Compressibility (1/Pa)
    dP = (pressure - P0) * 1e6  # Convert MPa to Pa
    
    # Calculate density
    rho = rho0 * (1 + beta * dP) / (1 + alpha * dT)
    
    print(f"DEBUG: Water density at T={temperature}K, P={pressure}MPa: {rho:.2f} kg/m³")
    return rho

def calculate_water_density_iapws(temperature, pressure):
    """
    IAPWS-IF97 density for T-array, P-array.
    Clamps P to [0.1, 100] MPa.
    Returns an array of densities.
    """
    # 1) Clamp P elementwise
    P_clamped = np.clip(pressure, 0.1, 100.0)

    # 2) Define a scalar helper
    def rho_scalar(T, P):
        return IAPWS97(T=T, P=P).rho

    # 3) Vectorize it
    vec_rho = np.vectorize(rho_scalar)

    # 4) Call over the grids
    return vec_rho(temperature, P_clamped)

def calculate_debye_huckel_slope(temperature, pressure, coeffs):
    """
    Calculate Debye-Hückel limiting law slope using Equation 4.
    
    Args:
        temperature (float): Temperature in K
        pressure (float): Pressure in MPa
        coeffs (dict): Dictionary with Debye-Hückel coefficients
        
    Returns:
        float: Debye-Hückel limiting law slope
    """
    delta_p = pressure - 0.1
    
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
    
    print(f"DEBUG: Debye-Hückel slope (A_v) at T={temperature}K, P={pressure}MPa: {A_v:.8e}")
    return A_v

def calculate_phi_v0(temperature, pressure, coeffs):
    """
    Calculate infinite dilution molal volume using coefficients.
    
    Args:
        temperature (float): Temperature in K
        pressure (float): Pressure in MPa
        coeffs (dict): Dictionary with phi_v0 coefficients
        
    Returns:
        float: Infinite dilution molal volume in cm³/mol
    """
    delta_p = pressure - 0.1  # Pressure difference from atmospheric
    
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
    
    # Combine all terms
    phi_v0 = term1 + term2 + term3
    
    print(f"DEBUG: ϕ⁰_v at T={temperature}K, P={pressure}MPa: {phi_v0:.6e} cm³/mol")
    return phi_v0

def calculate_interaction_parameter(temperature, pressure, ionic_strength, coeffs):
    """
    Calculate interaction parameter BMX.
    
    Args:
        temperature (float): Temperature in K
        pressure (float): Pressure in MPa
        ionic_strength (float): Ionic strength
        coeffs (dict): Dictionary with interaction coefficients
        
    Returns:
        float: Interaction parameter
    """
    delta_p = pressure - 0.1
    
    # Calculate B0 (temperature and pressure dependent)
    B0 = (coeffs.get('B00', 0) + 
         coeffs.get('B10', 0) * temperature + 
         coeffs.get('B20', 0) * temperature**2)
    
    if delta_p > 0:
        B0 += ((coeffs.get('B01', 0) + 
               coeffs.get('B11', 0) * temperature + 
               coeffs.get('B21', 0) * temperature**2) * delta_p)
        
        if 'B02' in coeffs:
            B0 += ((coeffs.get('B02', 0) + 
                   coeffs.get('B12', 0) * temperature + 
                   coeffs.get('B22', 0) * temperature**2) * delta_p**2)
    
    # Calculate B1 (concentration dependent term)
    B1 = coeffs.get('B1', 0)
    
    # Calculate BMX (Equation 5-6 from paper)
    B_MX = B0 + B1 * ionic_strength
    
    print(f"DEBUG: Interaction parameter B_MX at T={temperature}K, P={pressure}MPa, I={ionic_strength:.3f}: {B_MX:.8e}")
    return B_MX