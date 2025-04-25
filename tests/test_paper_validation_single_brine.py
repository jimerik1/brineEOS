import requests
import json
import sys

def test_appendix_a_example(api_url="http://localhost:5099/api/v1/calculate_density"):
    """
    Test the brine density API against the example in Appendix A of the paper.
    """
    print("Testing brine density API against Appendix A example...")
    
    # Parameters from Appendix A
    payload = {
        "brine_type": "mixed",
        "salt_composition": {
            "CaCl2": 24.0,
            "CaBr2": 25.3
        },
        "pressure_interval": [0.1, 0.2],  # Just one pressure point at 0.1 MPa
        "pressure_resolution": 0.1,
        "temperature_interval": [297, 298],  # Just one temperature point at 297 K
        "temperature_resolution": 1
    }
    
    try:
        # Send request to the API
        print(f"Sending request to {api_url}...")
        response = requests.post(api_url, json=payload)
        response.raise_for_status()
        
        # Parse the response
        data = response.json()
        
        # Find the density value in the response (handling possible key format differences)
        pressure_key = next((k for k in data["densities"].keys() 
                           if float(k) >= 0.09 and float(k) <= 0.11), None)
        if not pressure_key:
            print(f"ERROR: Could not find pressure key around '0.1' in response")
            print(f"Available pressure keys: {list(data['densities'].keys())}")
            return False
        
        temp_key = next((k for k in data["densities"][pressure_key].keys() 
                        if float(k) >= 296.5 and float(k) <= 297.5), None)
        if not temp_key:
            print(f"ERROR: Could not find temperature key around '297' in response")
            print(f"Available temperature keys: {list(data['densities'][pressure_key].keys())}")
            return False
        
        calculated_density = data["densities"][pressure_key][temp_key]
        expected_density = 1574.0  # From Appendix A
        
        # Calculate the difference
        absolute_diff = abs(calculated_density - expected_density)
        percent_diff = (absolute_diff / expected_density) * 100
        
        # Print the results
        print("\nValidation Test: Appendix A Example (24.0% CaCl2, 25.3% CaBr2 at 297K, 0.1MPa)")
        print(f"Expected density: {expected_density:.2f} kg/m³")
        print(f"Calculated density: {calculated_density:.2f} kg/m³")
        print(f"Absolute difference: {absolute_diff:.2f} kg/m³")
        print(f"Percent difference: {percent_diff:.2f}%")
        
        # Check if within reasonable tolerance (1%)
        if percent_diff <= 1.0:
            print("VALIDATION PASSED: Within 1% of expected value")
            return True
        else:
            print("VALIDATION FAILED: Difference exceeds 1%")
            return False
        
    except Exception as e:
        print(f"ERROR: {e}")
        if 'data' in locals():
            print(f"Response content: {json.dumps(data, indent=2)}")
        return False

# Also add a test for a single salt brine
def test_single_salt_brine(api_url="http://localhost:5000/api/v1/calculate_density"):
    """
    Test a single salt brine case (CaCl2 at reference conditions)
    """
    print("\nTesting single salt brine calculation...")
    
    # We'll use a common concentration of CaCl2 brine
    payload = {
        "brine_type": "CaCl2",
        "base_density": 1300.0,  # A standard reference density
        "pressure_interval": [0.1, 0.2],
        "pressure_resolution": 0.1,
        "temperature_interval": [298.15, 299.15],  # 25°C
        "temperature_resolution": 1
    }
    
    try:
        response = requests.post(api_url, json=payload)
        response.raise_for_status()
        data = response.json()
        
        # Extract the calculated density
        pressure_key = next(iter(data["densities"]))
        temp_key = next(iter(data["densities"][pressure_key]))
        calculated_density = data["densities"][pressure_key][temp_key]
        
        print(f"\nSingle Salt Test (CaCl2, base_density: 1300.0 kg/m³, at 298.15K, 0.1MPa)")
        print(f"Calculated density: {calculated_density:.2f} kg/m³")
        
        # For a single salt at reference conditions, the density should be close to base_density
        absolute_diff = abs(calculated_density - 1300.0)
        percent_diff = (absolute_diff / 1300.0) * 100
        print(f"Difference from base density: {absolute_diff:.2f} kg/m³ ({percent_diff:.2f}%)")
        
        if percent_diff <= 2.0:
            print("VALIDATION PASSED: Calculated density is within 2% of base density")
            return True
        else:
            print("VALIDATION FAILED: Calculated density differs significantly from base density")
            return False
            
    except Exception as e:
        print(f"ERROR: {e}")
        return False

if __name__ == "__main__":
    # Allow for custom API URL from command line
    api_url = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:5099/api/v1/calculate_density"
    
    # Run both tests
    test_appendix_a_example(api_url)
    test_single_salt_brine(api_url)