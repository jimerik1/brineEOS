# utils/validators.py
"""
Validation functions for API input
"""

def validate_input(data):
    """
    Validate API input data.
    
    Args:
        data (dict): Input data to validate
        
    Returns:
        list: List of validation errors, empty if valid
    """
    errors = []
    
    # Check required fields
    required_fields = ['brine_type']
    
    # If brine_type is 'mixed', salt_composition is required
    if 'brine_type' in data and data['brine_type'].lower() == 'mixed':
        required_fields.append('salt_composition')
    
    # For calculation, we need these parameters
    required_fields.extend([
        'pressure_interval',
        'pressure_resolution',
        'temperature_interval',
        'temperature_resolution'
    ])
    
    # If brine_type is 'single', we need base_density or salt_composition (to get salt type)
    if 'brine_type' in data and data['brine_type'].lower() == 'single':
        if 'base_density' not in data and 'salt_composition' not in data:
            errors.append("For 'single' brine type, either 'base_density' or 'salt_composition' must be provided")
    
    # Check required fields
    for field in required_fields:
        if field not in data:
            errors.append(f"Missing required field: {field}")
    
    if errors:
        return errors
    
    # Validate brine_type
    valid_brine_types = ['single', 'mixed', 'NaCl', 'KCl', 'CaCl2', 'CaBr2', 'ZnBr2', 'ZnCl2']
    if data['brine_type'] not in valid_brine_types:
        errors.append(f"Invalid brine_type: {data['brine_type']}. Valid types are: {', '.join(valid_brine_types)}")
    
    # Validate salt_composition for mixed brines
    if 'salt_composition' in data:
        composition = data['salt_composition']
        
        if not isinstance(composition, dict):
            errors.append("salt_composition must be a dictionary mapping salt names to weight percentages")
        else:
            valid_salts = ['NaCl', 'KCl', 'CaCl2', 'CaBr2', 'ZnBr2', 'ZnCl2']
            
            # Check for valid salt types
            for salt in composition:
                if salt not in valid_salts:
                    errors.append(f"Invalid salt type: {salt}. Valid types are: {', '.join(valid_salts)}")
            
            # Check weight percentages
            for salt, pct in composition.items():
                if not isinstance(pct, (int, float)):
                    errors.append(f"Weight percentage for {salt} must be a number")
                elif pct < 0:
                    errors.append(f"Weight percentage for {salt} cannot be negative")
            
            # Check total composition (if all are valid numbers)
            if all(isinstance(pct, (int, float)) for salt, pct in composition.items()):
                total_pct = sum(pct for salt, pct in composition.items() if pct > 0)
                if total_pct > 100:
                    errors.append(f"Total salt composition exceeds 100%: {total_pct}%")
    
    # Validate base_density if provided
    if 'base_density' in data:
        if not isinstance(data['base_density'], (int, float)):
            errors.append("base_density must be a number")
        elif data['base_density'] <= 0:
            errors.append("base_density must be positive")
    
    # Validate intervals
    if not isinstance(data['pressure_interval'], list) or len(data['pressure_interval']) != 2:
        errors.append("pressure_interval must be a list with two elements [min, max]")
    elif not all(isinstance(p, (int, float)) for p in data['pressure_interval']):
        errors.append("pressure_interval values must be numbers")
    elif data['pressure_interval'][0] >= data['pressure_interval'][1]:
        errors.append("pressure_interval: min must be less than max")
    elif data['pressure_interval'][0] < 0.1:
        errors.append("pressure_interval: min cannot be less than 0.1 MPa")
    elif data['pressure_interval'][1] > 150:
        errors.append("pressure_interval: max cannot exceed 150 MPa (model validity limit)")
    
    if not isinstance(data['pressure_resolution'], (int, float)):
        errors.append("pressure_resolution must be a number")
    elif data['pressure_resolution'] <= 0:
        errors.append("pressure_resolution must be positive")
    
    if not isinstance(data['temperature_interval'], list) or len(data['temperature_interval']) != 2:
        errors.append("temperature_interval must be a list with two elements [min, max]")
    elif not all(isinstance(t, (int, float)) for t in data['temperature_interval']):
        errors.append("temperature_interval values must be numbers")
    elif data['temperature_interval'][0] >= data['temperature_interval'][1]:
        errors.append("temperature_interval: min must be less than max")
    elif data['temperature_interval'][0] < 273.15:
        errors.append("temperature_interval: min cannot be less than 273.15 K (0°C)")
    elif data['temperature_interval'][1] > 447:
        errors.append("temperature_interval: max cannot exceed 447 K (174°C, model validity limit)")
    
    if not isinstance(data['temperature_resolution'], (int, float)):
        errors.append("temperature_resolution must be a number")
    elif data['temperature_resolution'] <= 0:
        errors.append("temperature_resolution must be positive")
    
    # Check for too many points (performance consideration)
    if len(errors) == 0:
        pressure_points = int((data['pressure_interval'][1] - data['pressure_interval'][0]) / data['pressure_resolution']) + 1
        temp_points = int((data['temperature_interval'][1] - data['temperature_interval'][0]) / data['temperature_resolution']) + 1
        
        total_points = pressure_points * temp_points
        if total_points > 10000:
            errors.append(f"Too many calculation points requested ({total_points}). "
                          f"Please use a larger resolution or smaller intervals.")
    
    return errors