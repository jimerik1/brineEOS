import json
import os

def convert_water_density_to_interpolation_json(filepath, output_filename="water_density.json"):
    """
    Reads water density data from a file, converts it to a JSON object
    in the format expected by the interpolation script, and saves it
    to a JSON file.

    Args:
        filepath (str): The path to the input text file.
        output_filename (str): The name of the output JSON file.
                                Defaults to "water_density.json".
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Find the start and end of each section
        pressure_start = -1
        temp_start = -1
        density_start = -1

        for i, line in enumerate(lines):
            if "PRESSURE VECTOR (Pa)" in line:
                pressure_start = i
            elif "TEMPERATURE VECTOR (Deg C)" in line:
                temp_start = i
            elif "WATER DENSITY (kg/m3)" in line:
                density_start = i

        if pressure_start == -1 or temp_start == -1 or density_start == -1:
            print("Error: Could not find all required sections (Pressure, Temperature, Density) in the file.")
            return

        # Extract and parse pressure data
        pressure_lines = []
        for i in range(pressure_start + 1, temp_start):
            pressure_lines.append(lines[i].strip())
        pressure_str = ' '.join(pressure_lines).replace('\t', ' ')
        # Convert pressures to float for keys, even though they are integers
        pressures = [float(p) for p in pressure_str.split()]

        # Extract and parse temperature data
        temp_lines = []
        for i in range(temp_start + 1, density_start):
            temp_lines.append(lines[i].strip())
        temp_str = ' '.join(temp_lines).replace('\t', ' ')
        # Replace commas with periods for proper float conversion
        temperatures = [float(t.replace(',', '.')) for t in temp_str.split()]

        # Extract and parse density data
        density_lines = []
        for i in range(density_start + 1, len(lines)):
            density_lines.append(lines[i].strip())
        density_str = ' '.join(density_lines).replace('\t', ' ')
        # Replace commas with periods for proper float conversion
        densities = [float(d.replace(',', '.')) for d in density_str.split()]

        # Organize the data into the desired JSON structure for interpolation
        json_data = {"densities": {}}
        temps_per_pressure = len(temperatures)

        if len(pressures) * temps_per_pressure != len(densities):
             print("Error: Mismatch between the number of pressure points, temperature points, and density values.")
             return

        density_index = 0
        for i, pressure in enumerate(pressures):
            pressure_key = str(pressure)
            json_data["densities"][pressure_key] = {}
            for j in range(temps_per_pressure):
                temp = temperatures[j]
                density = densities[density_index]
                temp_key = str(round(temp, 2)) # Use rounded temp as key
                json_data["densities"][pressure_key][temp_key] = density
                density_index += 1

        # Write the JSON data to a file
        with open(output_filename, 'w') as outfile:
            json.dump(json_data, outfile, indent=4)

        print(f"Successfully converted data and saved to '{output_filename}'")

    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Specify the path to your input file
file_path = 'waterdensity.txt'
output_file = 'water_density.json'

# Convert the data and save to a JSON file
convert_water_density_to_interpolation_json(file_path, output_file)
