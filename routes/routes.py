# routes.py
"""
API routes for brine density calculation service
"""
from flask import Blueprint, request, jsonify
from calculators import get_calculator
from utils.validators import validate_input

api_bp = Blueprint('api', __name__)

@api_bp.route('/api/v1/calculate_density', methods=['POST'])
def calculate_density():
    """
    Calculate brine density at various pressure and temperature conditions.
    
    Expected payload:
    {
        "brine_type": "single" | "mixed" | "NaCl" | "KCl" | "CaCl2" | "CaBr2" | "ZnBr2" | "ZnCl2",
        "salt_composition": {
            "NaCl": 0.0,     // Weight percentage
            "KCl": 0.0,      // Weight percentage
            "CaCl2": 0.0,    // Weight percentage
            "CaBr2": 0.0,    // Weight percentage
            "ZnBr2": 0.0,    // Weight percentage
            "ZnCl2": 0.0     // Weight percentage
        },
        "base_density": 1100.0,  // kg/m³ at reference conditions (298.15K, 0.1MPa)
        "pressure_interval": [0.1, 100.0],  // MPa
        "pressure_resolution": 10.0,  // MPa
        "temperature_interval": [298.15, 423.15],  // K
        "temperature_resolution": 10.0  // K
    }
    
    Returns:
    {
        "metadata": {
            "brine_type": string,
            "salt_composition": object,
            "base_density": number,
            "pressure_points": array,
            "temperature_points": array,
            "units": {
                "density": "kg/m³",
                "pressure": "MPa",
                "temperature": "K"
            }
        },
        "densities": {
            "pressure1": {
                "temp1": density,
                "temp2": density,
                ...
            },
            ...
        }
    }
    """
    try:
        data = request.get_json()
        
        # Validate input
        errors = validate_input(data)
        if errors:
            return jsonify({"error": "Validation error", "details": errors}), 400
        
        # Get appropriate calculator
        brine_type = data['brine_type']
        salt_composition = data.get('salt_composition')
        
        try:
            calculator = get_calculator(brine_type, salt_composition)
        except ValueError as e:
            return jsonify({'error': str(e)}), 400
        
        # Calculate densities
        try:
            if brine_type.lower() == 'mixed':
                results = calculator.calculate(
                    salt_composition=salt_composition,
                    pressure_interval=data['pressure_interval'],
                    pressure_resolution=data['pressure_resolution'],
                    temperature_interval=data['temperature_interval'],
                    temperature_resolution=data['temperature_resolution']
                )
            else:
                # Single salt case
                results = calculator.calculate(
                    base_density=data.get('base_density'),
                    pressure_interval=data['pressure_interval'],
                    pressure_resolution=data['pressure_resolution'],
                    temperature_interval=data['temperature_interval'],
                    temperature_resolution=data['temperature_resolution']
                )
            
            return jsonify(results)
            
        except Exception as e:
            return jsonify({'error': f'Calculation error: {str(e)}'}), 500
    
    except Exception as e:
        return jsonify({'error': f'Server error: {str(e)}'}), 500