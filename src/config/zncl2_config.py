"""
Configuration for ZnCl2 density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'ZnCl2',
    'molecular_weight': 0.13630,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Zn',
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
        'A00': 1.6470e-5,
        'A10': 0.0,  # No value in table
        'A20': 0.0,  # No value in table
        'A01': 0.0,  # No value in table
        'A11': 0.0,  # No value in table
        'A21': 0.0,  # No value in table
        'A02': 0.0,  # No value in table
        'A12': 0.0,  # No value in table
        'A22': 0.0   # No value in table
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 1.9141e-6,
        'B10': 0.0,  # No value in table
        'B20': 0.0,  # No value in table
        'B01': 0.0,  # No value in table
        'B11': 0.0,  # No value in table
        'B21': 0.0,  # No value in table
        'B02': 0.0,  # No value in table
        'B12': 0.0,  # No value in table
        'B1': -1.7875e-8
    }
}