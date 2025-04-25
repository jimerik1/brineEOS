# routes/water_routes.py
"""
API routes for pure water density calculation service
"""
from flask import Blueprint, request, jsonify
from calculators.water_calculator import WaterCalculator
from utils.validators import validate_water_input

water_api_bp = Blueprint('water_api', __name__, url_prefix='/api/v1')

@water_api_bp.route('/calculate_water_density', methods=['POST'])
def calculate_water_density_route():
    """
    Calculate pure water density at various pressure and temperature conditions using CoolProp.

    Expected payload:
    {
        "pressure_interval": [0.1, 100.0],  // MPa
        "pressure_resolution": 10.0,  // MPa
        "temperature_interval": [273.15, 423.15],  // K
        "temperature_resolution": 10.0  // K
    }

    Returns:
    {
        "metadata": {
            "fluid": "Water",
            "source": "CoolProp",
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
        errors = validate_water_input(data)
        if errors:
            return jsonify({"error": "Validation error", "details": errors}), 400

        # Get the calculator
        calculator = WaterCalculator()

        # Calculate densities
        try:
            results = calculator.calculate(
                pressure_interval=data['pressure_interval'],
                pressure_resolution=data['pressure_resolution'],
                temperature_interval=data['temperature_interval'],
                temperature_resolution=data['temperature_resolution']
            )

            return jsonify(results)

        except Exception as e:
            # Log the error if needed
            # current_app.logger.error(f'Water calculation error: {e}')
            return jsonify({'error': f'Calculation error: {str(e)}'}), 500

    except Exception as e:
        # Log the error if needed
        # current_app.logger.error(f'Water endpoint server error: {e}')
        return jsonify({'error': f'Server error: {str(e)}'}), 500
