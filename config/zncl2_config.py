# config/znclc_config.py
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
        'A00': 1.1470e-5,  # Limited data available
        'A10': 0.0,
        'A20': 0.0,
        'A01': 0.0,
        'A11': 0.0,
        'A21': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 1.8141e-6,
        'B10': 0.0,
        'B20': 0.0,
        'B01': 0.0,
        'B11': 0.0,
        'B21': 0.0,
        'B1': -1.7876e-8
    }
}