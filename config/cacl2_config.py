# config/cacl2_config.py
"""
Configuration for CaCl2 density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'CaCl2',
    'molecular_weight': 0.11099,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Ca',
            'charge': 2,
            'stoichiometry': 1
        },
        'anion': {
            'symbol': 'Cl',
            'charge': -1,
            'stoichiometry': 2
        }
    },
    'phi_v0_coeffs': {
        'A00': -13.4695e-5,
        'A10': 9.824e-7,
        'A20': -15.783e-10,
        'A01': 1.348e-7,
        'A11': 48.2362e-10,
        'A21': -27.4851e-12,
        'A02': 34.3983e-10,
        'A12': 0.0,
        'A22': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 4.9949e-6,
        'B10': -22.6474e-9,
        'B20': 3.3664e-11,
        'B01': -7.1662e-8,
        'B11': 2.8029e-11,
        'B21': 0.0,
        'B02': 0.0,
        'B12': -2.2953e-13,
        'B22': 0.0,
        'B1': -2.2953e-8
    }
}