"""
Configuration for ZnBr2 density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'ZnBr2',
    'molecular_weight': 0.22536,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Zn',
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
        'A00': -11.0575e5,
        'A10': 7.5860e7,
        'A20': -10.1729e10,
        'A01': 0.0, # No value in table
        'A11': -3.0987e10,
        'A21': 0.0,  # No value in table
        'A02': 0.0,  # No value in table
        'A12': 0.0,  # No value in table
        'A22': 5.8225e15
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 5.5069e6,
        'B10': -6.1731e9,
        'B20': 0.0,  # No value in table
        'B01': 0.0,  # No value in table
        'B11': 2.7980e11,
        'B21': 0.0,
        'B02': 0.0,  # No value in table
        'B12': -1.7568e13,
        'B1': -9.2343e8,
    }
}