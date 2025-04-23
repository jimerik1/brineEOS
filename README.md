Brine Density API
A Flask-based REST API for calculating brine densities at various pressures and temperatures, based on the methodology described in SPE/IADC 16079: "Density Modeling for Pure and Mixed-Salt Brines as a Function of Composition, Temperature, and Pressure."
Table of Contents

Overview
Features
Installation
Usage
API Reference
Docker Deployment
Project Structure
Background
License

Overview
This API implements the density model from SPE/IADC 16079 to calculate densities for various brine types across a range of temperatures and pressures. It supports both single-salt brines (NaCl, KCl, CaCl₂, CaBr₂, ZnBr₂, ZnCl₂) and mixed-salt brines.
The calculations account for:

Salt composition effects
Temperature dependence (273.15K to 447K)
Pressure effects (0.1MPa to 150MPa)
Ion interactions in mixed brines

Features

Calculates density for single-salt brines using base density or composition
Handles complex mixed-salt brine calculations
Comprehensive input validation
JSON-based API with detailed metadata
Docker support for easy deployment
Health check endpoint for monitoring

Installation
Prerequisites

Python 3.9+
pip (Python package manager)

Local Installation

Clone the repository:
bashgit clone https://github.com/yourusername/brine-density-api.git
cd brine-density-api

Create a virtual environment (recommended):
bashpython -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

Install dependencies:
bashpip install -r requirements.txt

Run the API:
bashpython app.py


The API will be available at http://localhost:5000.
Usage
Example: Calculate NaCl Brine Density
pythonimport requests
import json

url = "http://localhost:5000/api/v1/calculate_density"

payload = {
    "brine_type": "NaCl",
    "base_density": 1100.0,
    "pressure_interval": [0.1, 100.0],
    "pressure_resolution": 25.0,
    "temperature_interval": [298.15, 373.15],
    "temperature_resolution": 25.0
}

response = requests.post(url, json=payload)
data = response.json()

print(json.dumps(data, indent=2))
Example: Calculate Mixed Brine Density
pythonimport requests
import json

url = "http://localhost:5000/api/v1/calculate_density"

payload = {
    "brine_type": "mixed",
    "salt_composition": {
        "CaCl2": 24.0,
        "CaBr2": 25.3
    },
    "pressure_interval": [0.1, 50.0],
    "pressure_resolution": 10.0,
    "temperature_interval": [298.15, 348.15],
    "temperature_resolution": 10.0
}

response = requests.post(url, json=payload)
data = response.json()

print(json.dumps(data, indent=2))
API Reference
Calculate Density
Endpoint: POST /api/v1/calculate_density
Calculate brine densities at various pressure and temperature conditions.
Request Body:
json{
  "brine_type": "single" | "mixed" | "NaCl" | "KCl" | "CaCl2" | "CaBr2" | "ZnBr2" | "ZnCl2",
  "salt_composition": {
    "NaCl": 0.0,
    "KCl": 0.0,
    "CaCl2": 0.0,
    "CaBr2": 0.0,
    "ZnBr2": 0.0,
    "ZnCl2": 0.0
  },
  "base_density": 1100.0,
  "pressure_interval": [0.1, 100.0],
  "pressure_resolution": 10.0,
  "temperature_interval": [298.15, 423.15],
  "temperature_resolution": 10.0
}
Parameters:
ParameterTypeRequiredDescriptionbrine_typestringYesType of brine: "single", "mixed", or a specific salt namesalt_compositionobjectYes for "mixed"Weight percentages of salts in the brinebase_densitynumberFor "single" if no compositionBase density at reference conditions (kg/m³)pressure_intervalarrayYesMin and max pressure in MPa [min, max]pressure_resolutionnumberYesPressure step size in MPatemperature_intervalarrayYesMin and max temperature in K [min, max]temperature_resolutionnumberYesTemperature step size in K
Response:
json{
  "metadata": {
    "brine_type": "NaCl",
    "base_density": 1100.0,
    "pressure_points": [0.1, 25.1, 50.1, 75.1, 100.1],
    "temperature_points": [298.15, 323.15, 348.15, 373.15],
    "units": {
      "density": "kg/m³",
      "pressure": "MPa",
      "temperature": "K"
    }
  },
  "densities": {
    "0.10": {
      "298.15": 1100.00,
      "323.15": 1090.25,
      "348.15": 1079.83,
      "373.15": 1068.76
    },
    "25.10": {
      "298.15": 1110.32,
      "323.15": 1100.68,
      "348.15": 1090.37,
      "373.15": 1079.42
    },
    // Additional pressure/temperature points...
  }
}
Health Check
Endpoint: GET /health
Check the status of the API.
Response:
json{
  "status": "ok"
}
Docker Deployment
Build and Run with Docker
bash# Build the Docker image
docker build -t brine-density-api .

# Run the container
docker run -d -p 5000:5000 --name brine-density-api brine-density-api
Using Docker Compose
bash# Start the service
docker-compose up -d

# View logs
docker-compose logs -f

# Stop the service
docker-compose down
Project Structure
project/
│
├── app.py                  # Main Flask application
├── routes.py               # API endpoints
├── calculators/            
│   ├── __init__.py         # Factory function
│   ├── single_brine_calculator.py  # Single salt calculator
│   └── mixed_brine_calculator.py   # Mixed salt calculator
│
├── config/                 # Configuration files with parameters from paper
│   ├── __init__.py
│   ├── constants.py        # Physical constants
│   ├── nacl_config.py      # Parameters for NaCl
│   ├── kcl_config.py       # Parameters for KCl
│   ├── cacl2_config.py     # Parameters for CaCl₂
│   └── ... (other salts)
│
├── utils/                  # Utility functions
│   ├── __init__.py
│   └── validators.py       # Input validators
│
├── requirements.txt        # Python dependencies
├── Dockerfile              # Docker build instructions
└── docker-compose.yml      # Docker Compose configuration
Background
This API implements the density model described in SPE/IADC 16079 paper "Density Modeling for Pure and Mixed-Salt Brines as a Function of Composition, Temperature, and Pressure" by N.P. Kemp and D.C. Thomas.
The model uses the Brønsted-Guggenheim extension of the Debye-Hückel theory to calculate apparent molal volumes of electrolytes in solution. It accounts for ion-ion interactions and their dependency on temperature and pressure.
Key aspects of the model:

Calculation of infinite dilution molal volumes (φ°v)
Debye-Hückel term for long-range electrostatic interactions
Ion-specific interaction parameters for short-range interactions
Normalization factors for mixed salt brines

