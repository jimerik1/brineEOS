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
        'A00': -13.4695e-5,
        'A10': 9.8240e-7,
        'A20': -15.7830e-10,
        'A01': 0.0, # No value in table
        'A11': 0.0, # No value in table
        'A21': 1.2489e-12,
        'A02': 48.2352e-10,
        'A12': -27.4851e-12,
        'A22': 34.3693e-15,

    },
    'debye_huckel_coeffs': DEBYE_HUCKEL_COEFFS,
    'interaction_coeffs': {
        'B00': 4.9949e-6,
        'B10': -22.6474e-9,
        'B20': 3.3864e-11,
        'B01': -7.1552e-9,
        'B11': 0.0,  # No value in table
        'B21': 0.0,  # No value in table
        'B02': 2.8026e-11,
        'B12': 0.0,  # No value in table
        'B1': -2.2953e-8
    }
}