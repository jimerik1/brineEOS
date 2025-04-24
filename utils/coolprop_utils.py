# utils/coolprop_utils.py
from CoolProp.CoolProp import PropsSI

def calculate_water_density_cp(temperature: float, pressure: float, phase: str = None) -> float:
    """
    Return pure-water density [kg/m³] via CoolProp.
    
    Args:
        temperature: in K
        pressure: in MPa
        phase: None (default) → returns whichever phase CoolProp gives,
               'liquid' → saturated liquid at that P,
               'vapor'  → saturated vapor at that P,
               'tp'     → two-phase interpolated by T/P.
    """
    P_Pa = pressure * 1e6
    if phase == 'liquid':
        # saturated liquid at P
        return PropsSI('Dmass', 'Q', 0, 'P', P_Pa, 'Water')
    elif phase == 'vapor':
        # saturated vapor at P
        return PropsSI('Dmass', 'Q', 1, 'P', P_Pa, 'Water')
    else:
        # default: pure T–P lookup (may jump to vapor if T exceeds sat. temp)
        return PropsSI('Dmass', 'T', temperature, 'P', P_Pa, 'Water')