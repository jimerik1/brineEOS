"""
Configuration for KCl density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'KCl',
    'molecular_weight': 0.07455,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'K',
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
        'A00': -6.0780e5,
        'A10': 5.4193e7,
        'A20': -8.2782e10,
        'A01': 0.4514e7,
        'A11': 0.0,  # No value in table
        'A21': 0.0,  # No value in table
        'A02': -1.5741e10,
        'A12': 0.0,  # No value in table
        'A22': 0.0   # No value in table
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 9.1287e6,
        'B10': -49.3068e9,
        'B20': 7.3418e11,
        'B01': -4.1250e8,
        'B11': 0.0,  # No value in table
        'B21': 0.0,  # No value in table
        'B02': 0.0,  # No value in table
        'B12': 0.0,  # No value in table
        'B1': -7.1118e8
    }
}