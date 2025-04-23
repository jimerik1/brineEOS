# config/cabr2_config.py (continued)
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
        'A00': -13.5882e-5,
        'A10': 10.5681e-7,
        'A20': -16.5140e-10,
        'A01': 1.1752e-7,
        'A11': 52.5124e-10,
        'A21': -28.5896e-12,
        'A02': 37.3889e-10,
        'A12': 0.0,
        'A22': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 6.8218e-6,
        'B10': -25.4804e-9,
        'B20': 0.0,
        'B01': -7.4830e-8,
        'B11': 2.2943e-11,
        'B21': 0.0,
        'B02': 0.0,
        'B12': -3.4634e-13,
        'B22': 0.0,
        'B1': 0.0
    }
}