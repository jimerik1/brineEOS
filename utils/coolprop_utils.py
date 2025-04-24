# utils/coolprop_utils.py
from CoolProp.CoolProp import PropsSI

def calculate_water_density_cp(temperature: float, pressure: float) -> float:
    """
    Return pure-water density [kg/m³] via CoolProp.

    Args:
        temperature: in K
        pressure: in MPa
    """
    # PropsSI expects pressure in Pa
    P_Pa = pressure * 1e6
    # 'D' = mass density; inputs: T (K), P (Pa), fluid name
    rho = PropsSI('D', 'T', temperature, 'P', P_Pa, 'Water')  #  [oai_citation_attribution:0‡coolprop.org](https://www.coolprop.org/coolprop/HighLevelAPI.html?utm_source=chatgpt.com)
    return rho