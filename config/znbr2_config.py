# config/znbr2_config.py
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
        'A00': -11.0575e-5,
        'A10': 7.5860e-7,
        'A20': -10.1728e-10,
        'A01': -3.0887e-7,
        'A11': 0.0,
        'A21': 0.0,
        'A02': 5.8225e-10,
        'A12': 0.0,
        'A22': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 5.6088e-6,
        'B10': -8.1731e-9,
        'B20': 0.0,
        'B01': 2.7880e-8,
        'B11': 0.0,
        'B21': 0.0,
        'B02': -1.7688e-11,
        'B12': -8.2343e-13,
        'B22': 0.0,
        'B1': 0.0
    }
}