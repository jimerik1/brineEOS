"""
Configuration for CaBr2 density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'CaBr2',
    'molecular_weight': 0.19989,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Ca',
            'charge': 2,
            'stoichiometry': 1
        },
        'anion': {
            'symbol': 'Br',
            'charge': -1,
            'stoichiometry': 2
        }
    },
    'phi_v0_coeffs': {
        'A00': -13.6982e5,
        'A10': 10.5661e7,
        'A20': -16.5140e10,
        'A01': 0.0,  # No value in table
        'A11': 0.0,  # No value in table
        'A21': 1.1752e12,
        'A02': 52.5124e10,
        'A12': -29.5866e12,
        'A22': 37.3999e15,
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 5.8218e6,
        'B10': -25.4804e9,
        'B20': 3.6141e11,
        'B01': -7.4830e9,
        'B11': 0.0,  # No value in table
        'B21': 0.0,  # No value in table
        'B02': 2.2843e11,
        'B12': 0.0,  # No value in table
        'B1': -3.4634e8
    }
}