# config/nacl_config.py
"""
Configuration for NaCl density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'NaCl',
    'molecular_weight': 0.05844,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Na',
            'charge': 1,
            'stoichiometry': 1
        },
        'anion': {
            'symbol': 'Cl',
            'charge': -1,
            'stoichiometry': 1
        }
    },
    'phi_v0_coeffs': {
        'A00': -7.1828e-5,
        'A10': 5.8004e-7,
        'A20': -8.8181e-10,
        'A01': 5.3789e-7,
        'A11': -32.1810e-10,
        'A21': 4.7183e-12,
        'A02': 2.3388e-10,
        'A12': -3.8432e-12,
        'A22': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 2.7167e-6,
        'B10': 0.0,
        'B20': 0.0,
        'B01': 0.0,
        'B11': 7.8163e-14,
        'B21': 0.0,
        'B02': 0.0,
        'B12': 3.6462e-13,
        'B22': 0.0,
        'B1': -9.0242e-8
    }
}