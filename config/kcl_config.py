# config/kcl_config.py
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
        'A00': -6.0780e-5,
        'A10': 6.4183e-7,
        'A20': -9.2782e-10,
        'A01': 4.4614e-7,
        'A11': -1.6741e-11,
        'A21': 0.0
    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 8.1287e-6,
        'B10': -48.3096e-9,
        'B20': 0.0,
        'B01': -4.1260e-8,
        'B11': 0.0,
        'B21': 0.0,
        'B1': -7.1118e-8
    }
}