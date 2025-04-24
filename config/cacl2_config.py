"""
Configuration for CaCl2 density calculations from Tables 2 and 3
"""
from config.constants import DEBYE_HUCKEL_COEFFS

CONFIG = {
    'name': 'CaCl2',
    'molecular_weight': 0.11099,  # kg/mol
    'ions': {
        'cation': {
            'symbol': 'Ca',
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
        'A00': -13.4695e5,
        'A10': 9.8240e7,
        'A20': -15.7830e10,
        'A01': 0.0, # No value in table
        'A11': 0.0, # No value in table
        'A21': 1.2489e12,
        'A02': 48.2352e10,
        'A12': -27.4851e12,
        'A22': 34.3693e15,

    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 4.9949e6,
        'B10': -22.6474e9,
        'B20': 3.3864e11,
        'B01': -7.1552e9,
        'B11': 0.0,  # No value in table
        'B21': 0.0,  # No value in table
        'B02': 2.8026e11,
        'B12': 0.0,  # No value in table
        'B1': -2.2953e8
    }
}