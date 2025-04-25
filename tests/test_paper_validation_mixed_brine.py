#!/usr/bin/env python3
"""
Brine API Validation Script

This script validates the Brine Density API by sending multiple payloads
and plotting the results as points without connecting lines.
"""

import requests
import json
import matplotlib.pyplot as plt
import numpy as np

# API endpoint - using port 5099 as confirmed working
API_URL = "http://localhost:5099/api/v1/calculate_density"

# Define the payloads
payloads = [
    {
        "brine_type": "mixed",
        "salt_composition": {
            "CaCl2": 24.0,
            "CaBr2": 25.3
        },
        "pressure_interval": [1.0, 150.0],
        "pressure_resolution": 1.0,
        "temperature_interval": [295.00, 296.00],
        "temperature_resolution": 1.0
    },
    {
        "brine_type": "mixed",
        "salt_composition": {
            "CaCl2": 24.0,
            "CaBr2": 25.3
        },
        "pressure_interval": [1.0, 150.0],
        "pressure_resolution": 1.0,
        "temperature_interval": [366.00, 367.00],
        "temperature_resolution": 1.0
    },
    {
        "brine_type": "mixed",
        "salt_composition": {
            "CaCl2": 24.0,
            "CaBr2": 25.3
        },
        "pressure_interval": [1.0, 150.0],
        "pressure_resolution": 1.0,
        "temperature_interval": [449.00, 450.00],
        "temperature_resolution": 1.0
    }
]

# Labels for the plot
labels = ["295K", "366K", "449K"]

# Check if the last payload might be invalid
if payloads[2]["temperature_interval"][0] > 447:
    print(f"Warning: Temperature {payloads[2]['temperature_interval'][0]}K exceeds model validity limit (447K)")
    print("Adjusting to valid range...")
    payloads[2]["temperature_interval"] = [446.00, 447.00]
    labels[2] = "446K"

def send_request(payload):
    """Send a request to the API and return the response"""
    try:
        print(f"Sending request to {API_URL}")
        response = requests.post(API_URL, json=payload)
        
        # Print status code for debugging
        print(f"Response status code: {response.status_code}")
        
        # If there's an error status code, print the response text for debugging
        if response.status_code >= 400:
            print(f"Error response: {response.text}")
            
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error sending request: {e}")
        return None

def extract_densities(response):
    """Extract pressure and density values from the API response"""
    if not response or "densities" not in response:
        print("No 'densities' field in response or response is None")
        if response:
            print(f"Response keys: {response.keys()}")
        return [], []
    
    pressures = []
    densities = []
    
    # Sort pressure points to ensure they're in ascending order
    sorted_pressure_keys = sorted([float(p) for p in response["densities"].keys()])
    
    # Extract the first temperature key from each pressure point
    for pressure in sorted_pressure_keys:
        pressure_str = f"{pressure:.2f}"
        temp_data = response["densities"].get(pressure_str)
        
        if temp_data:
            # Get the first temperature key
            temp_keys = list(temp_data.keys())
            if temp_keys:
                temp_key = temp_keys[0]
                density = temp_data[temp_key]
                
                pressures.append(pressure)
                densities.append(density)
    
    return pressures, densities

def plot_results(all_data):
    """Plot the results as points only (no connecting lines)"""
    plt.figure(figsize=(12, 8))
    
    # Collect all pressures and densities into flat lists
    all_pressures = []
    all_densities = []
    for pressures, densities, _ in all_data:
        all_pressures.extend(pressures)
        all_densities.extend(densities)
    
    # Plot each series
    markers = ['o', 's', '^']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    for i, (pressures, densities, label) in enumerate(all_data):
        if pressures and densities:
            plt.scatter(
                pressures, densities,
                marker=markers[i % len(markers)],
                color=colors[i % len(colors)],
                s=20,
                alpha=0.7,
                label=label
            )
    
    # Dynamic axis limits with 5% padding
    x_min, x_max = min(all_pressures), max(all_pressures)
    y_min, y_max = min(all_densities), max(all_densities)
    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05
    plt.xlim(x_min - x_pad, x_max + x_pad)
    plt.ylim(y_min - y_pad, y_max + y_pad)
    
    plt.title('Brine Density vs Pressure for CaCl2(24%) + CaBr2(25.3%) Mixture')
    plt.xlabel('Pressure (MPa)')
    plt.ylabel('Density (kg/m³)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('brine_density_points_only.png', dpi=300)
    print("Plot saved as 'brine_density_points_only.png'")
    plt.show()
    
def main():
    """Main function to run the validation"""
    all_data = []
    
    # First check if API is accessible
    try:
        health_url = API_URL.rsplit('/', 2)[0] + '/health'
        health_response = requests.get(health_url)
        if health_response.status_code == 200:
            print(f"API health check successful: {health_response.json()}")
        else:
            print(f"API health check failed with status code: {health_response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"API health check failed: {e}")
        print("Make sure the API is running and accessible at the correct port (5099)")
        return
    
    for i, payload in enumerate(payloads):
        print(f"\nProcessing payload {i+1}:")
        print(f"Temperature interval: {payload['temperature_interval']}K")
        
        response = send_request(payload)
        if not response:
            print(f"No response received for payload {i+1}")
            all_data.append(([], [], labels[i]))
            continue
        
        pressures, densities = extract_densities(response)
        if not pressures or not densities:
            print(f"No valid data extracted for payload {i+1}")
            all_data.append(([], [], labels[i]))
            continue
            
        print(f"Data points: {len(pressures)}")
        print(f"Density range: {min(densities):.2f} - {max(densities):.2f} kg/m³")
        
        all_data.append((pressures, densities, labels[i]))
    
    if any(len(pressures) > 0 for pressures, _, _ in all_data):
        plot_results(all_data)
    else:
        print("\nNo valid data to plot. Please check the API connection and payload validity.")

if __name__ == "__main__":
    main()