# utils/conversions.py
def to_mpa(value: float, unit: str) -> float:
    if unit.lower() == "mpa":
        return value
    if unit.lower() == "bar":
        return value / 10.0          # 1 MPa = 10 bar
    raise ValueError(f"Unsupported pressure unit '{unit}'")

def to_kelvin(value: float, unit: str) -> float:
    if unit.lower() == "k":
        return value
    if unit.lower() == "c":
        return value + 273.15
    raise ValueError(f"Unsupported temperature unit '{unit}'")