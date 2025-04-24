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
        'A00': -7.1829e5,
        'A10': 5.6004e7,
        'A20': -8.6619e10,
        'A01': 5.3789e7,
        'A11': -32.1990e10,
        'A21': 4.7893e12,
        'A02': 0.0,  # No value in table
        'A12': 2.3368e12,
        'A22': -3.9432e15
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 0.0,  # No value in table
        'B10': 2.7157e9,
        'B20': 0.0,  # No value in table
        'B01': 0.0,  # No value in table 
        'B11': 0.0,  # No value in table
        'B21': 7.6153e14,
        'B02': 0.0,  # No value in table
        'B12': 3.5482e13,
        'B1': -9.0242e8
    }
}