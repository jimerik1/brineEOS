import requests
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_atmospheric_water_table(table_text):
    """Parse the atmospheric pressure water density table text into temperature and density arrays."""
    data = []
    
    lines = table_text.strip().split('\n')
    
    # Find the start of the data table
    start_idx = 0
    for i, line in enumerate(lines):
        if '°C' in line and '.0' in line and '.1' in line:
            start_idx = i + 1
            break
    
    # Parse the data
    for line in lines[start_idx:]:
        if not line.strip() or 'Literature sources' in line:
            break
        
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        
        try:
            temp_base = float(parts[0])
            for i in range(1, len(parts)):
                if parts[i].strip():
                    temp_c = temp_base + (i-1)/10
                    density = float(parts[i])
                    data.append((temp_c, density))
        except (ValueError, IndexError) as e:
            logging.warning(f"Error parsing line: {line}, error: {e}")
    
    # Convert to numpy arrays
    temps_c = np.array([d[0] for d in data])
    densities = np.array([d[1] for d in data])
    
    return temps_c, densities

def parse_pressure_water_table(table_text):
    """Parse the medium/high pressure water density table text into temperature and density arrays."""
    data = []
    
    lines = table_text.strip().split('\n')
    
    # Skip the header row
    for line in lines[1:]:  # Skip header row
        if not line.strip() or "boiling point" in line:  # Skip empty lines and the boiling point row
            continue
            
        try:
            parts = line.strip().split('\t')
            if len(parts) < 4:  # We need at least temperature in C and density in kg/m3
                continue
                
            temp_c = float(parts[0])
            density_kg_m3 = float(parts[3])  # Density in kg/m3 is in column 4
            
            data.append((temp_c, density_kg_m3))
                
        except (ValueError, IndexError) as e:
            logging.warning(f"Error parsing line: {line}, error: {e}")
    
    # Convert to numpy arrays
    temps_c = np.array([d[0] for d in data])
    densities = np.array([d[1] for d in data])
    
    return temps_c, densities

def test_water_api(temps_c, ref_densities, pressure_mpa, max_temp_c=100.0):
    """Test water density API endpoint against reference values at specified pressure."""
    api_url = "http://localhost:5099/api/v1/calculate_water_density"
    
    # Convert Celsius to Kelvin (temperatures in the API need to be in Kelvin)
    temps_k = temps_c + 273.15
    
    # Pressure name for logging
    if pressure_mpa < 1:
        pressure_name = "atmospheric"
    elif pressure_mpa < 10:
        pressure_name = "medium"
    else:
        pressure_name = "high"
    
    # Prepare arrays for results
    valid_temps_c = []
    api_densities = []
    valid_ref_densities = []
    
    # Test each temperature point individually
    for i, temp_k in enumerate(temps_k):
        # Skip temperatures higher than max_temp_c
        if temps_c[i] > max_temp_c:
            continue
            
        try:
            payload = {
                "pressure_interval": [pressure_mpa, pressure_mpa + 0.0001],  # Small range to ensure valid interval
                "pressure_resolution": 0.0001,
                "temperature_interval": [temp_k, temp_k + 0.0001],  # Small range to ensure valid interval
                "temperature_resolution": 0.0001
            }
            
            logging.info(f"Testing temp: {temps_c[i]:.1f}°C ({temp_k:.2f}K) at {pressure_mpa:.6f} MPa ({pressure_name} pressure)")
            response = requests.post(api_url, json=payload)
            
            if response.status_code != 200:
                logging.warning(f"Error for T={temp_k:.2f}K ({temps_c[i]:.1f}°C): {response.status_code} - {response.text}")
                continue
                
            data = response.json()
            
            # Find the closest pressure and temperature keys
            p_keys = list(data["densities"].keys())
            if not p_keys:
                logging.warning(f"No pressure data found for T={temp_k:.2f}K")
                continue
                
            p_key = p_keys[0]  # Just use the first pressure point
            t_keys = list(data["densities"][p_key].keys())
            if not t_keys:
                logging.warning(f"No temperature data found for P={p_key}")
                continue
                
            t_key = t_keys[0]  # Just use the first temperature point
            
            # Get the density
            density = data["densities"][p_key][t_key]
            if density is None:
                logging.warning(f"Null density value for T={t_key}K, P={p_key}MPa")
                continue
                
            # If we got here, we have a valid data point
            valid_temps_c.append(temps_c[i])
            api_densities.append(density)
            valid_ref_densities.append(ref_densities[i])
            
            # For atmospheric pressure, convert from g/ml to kg/m³ for logging
            if pressure_name == "atmospheric":
                ref_value = ref_densities[i] * 1000
            else:
                ref_value = ref_densities[i]
                
            logging.info(f"  Ref: {ref_value:.2f} kg/m³, API: {density:.2f} kg/m³")
            
        except Exception as e:
            logging.error(f"Exception for T={temp_k:.2f}K: {e}")
    
    # Convert to numpy arrays
    valid_temps_c = np.array(valid_temps_c)
    api_densities = np.array(api_densities)
    valid_ref_densities = np.array(valid_ref_densities)
    
    # For atmospheric pressure, need to convert from g/ml to kg/m³
    if pressure_name == "atmospheric":
        valid_ref_densities_kg_m3 = valid_ref_densities * 1000
    else:
        valid_ref_densities_kg_m3 = valid_ref_densities
    
    # Calculate differences
    diff = api_densities - valid_ref_densities_kg_m3
    
    # Only proceed if we have valid data
    if len(valid_temps_c) == 0:
        logging.error(f"No valid data points collected for {pressure_name} pressure.")
        return None
    
    # Calculate statistics
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    max_diff = np.max(np.abs(diff))
    
    print(f"\n{pressure_name.capitalize()} pressure water density comparison results:")
    print(f"Number of valid data points: {len(valid_temps_c)}")
    print(f"Mean difference: {mean_diff:.2f} kg/m³")
    print(f"Standard deviation: {std_diff:.2f} kg/m³")
    print(f"Maximum absolute difference: {max_diff:.2f} kg/m³")
    
    return {
        'temps_c': valid_temps_c,
        'ref_densities': valid_ref_densities_kg_m3,
        'api_densities': api_densities,
        'diff': diff,
        'pressure': pressure_mpa,
        'pressure_name': pressure_name
    }

def main():
    # Get the atmospheric water density table reference data
    atmospheric_water_table = """
°C	.0	.1	.2	.3	.4	.5	.6	.7	.8	.9
0	0.999843	0.999850	0.999856	0.999863	0.999869	0.999875	0.999880	0.999886	0.999892	0.999897
1	0.999902	0.999907	0.999911	0.999916	0.999920	0.999925	0.999929	0.999932	0.999936	0.999940
2	0.999943	0.999946	0.999949	0.999952	0.999955	0.999957	0.999960	0.999962	0.999964	0.999965
3	0.999967	0.999969	0.999970	0.999971	0.999972	0.999973	0.999974	0.999974	0.999975	0.999975
4	0.999975	0.999975	0.999974	0.999974	0.999973	0.999973	0.999972	0.999971	0.999970	0.999968
5	0.999967	0.999965	0.999963	0.999961	0.999959	0.999957	0.999954	0.999952	0.999949	0.999946
6	0.999943	0.999940	0.999936	0.999933	0.999929	0.999925	0.999922	0.999917	0.999913	0.999909
7	0.999904	0.999900	0.999895	0.999890	0.999885	0.999879	0.999874	0.999868	0.999863	0.999857
8	0.999851	0.999845	0.999839	0.999832	0.999826	0.999819	0.999812	0.999805	0.999798	0.999791
9	0.999784	0.999776	0.999768	0.999761	0.999753	0.999745	0.999737	0.999728	0.999720	0.999711
10	0.999702	0.999694	0.999685	0.999675	0.999666	0.999657	0.999647	0.999638	0.999628	0.999618
11	0.999608	0.999598	0.999587	0.999577	0.999566	0.999556	0.999545	0.999534	0.999523	0.999512
12	0.999500	0.999489	0.999477	0.999466	0.999454	0.999442	0.999430	0.999417	0.999405	0.999393
13	0.999380	0.999367	0.999355	0.999342	0.999328	0.999315	0.999302	0.999288	0.999275	0.999261
14	0.999247	0.999233	0.999219	0.999205	0.999191	0.999176	0.999162	0.999147	0.999133	0.999118
15	0.999103	0.999087	0.999072	0.999057	0.999041	0.999026	0.999010	0.998994	0.998978	0.998962
16	0.998946	0.998930	0.998913	0.998897	0.998880	0.998863	0.998847	0.998830	0.998813	0.998795
17	0.998778	0.998761	0.998743	0.998725	0.998708	0.998690	0.998672	0.998654	0.998635	0.998617
18	0.998599	0.998580	0.998561	0.998543	0.998524	0.998505	0.998486	0.998467	0.998447	0.998428
19	0.998408	0.998389	0.998369	0.998349	0.998329	0.998309	0.998289	0.998269	0.998248	0.998228
20	0.998207	0.998186	0.998166	0.998145	0.998124	0.998103	0.998081	0.998060	0.998039	0.998017
25	0.997048	0.997022	0.996996	0.996970	0.996944	0.996918	0.996892	0.996866	0.996839	0.996813
30	0.99565	0.99562	0.99559	0.99556	0.99553	0.99550	0.99547	0.99544	0.99541	0.99537
35	0.99403	0.99400	0.99396	0.99393	0.99390	0.99386	0.99383	0.99379	0.99376	0.99372
40	0.99222	0.99218	0.99214	0.99210	0.99206	0.99202	0.99199	0.99195	0.99191	0.99187
45	0.99021	0.99017	0.99013	0.99009	0.99004	0.99000	0.98996	0.98992	0.98988	0.98983
50	0.98804	0.98799	0.98794	0.98790	0.98785	0.98781	0.98776	0.98772	0.98767	0.98763
55	0.98569	0.98564	0.98560	0.98555	0.98550	0.98545	0.98540	0.98535	0.98530	0.98525
60	0.98320	0.98314	0.98309	0.98304	0.98299	0.98294	0.98289	0.98283	0.98278	0.98273
65	0.98055	0.98050	0.98044	0.98039	0.98033	0.98028	0.98022	0.98017	0.98011	0.98006
70	0.97776	0.97771	0.97765	0.97759	0.97754	0.97748	0.97742	0.97736	0.97731	0.97725
75	0.97484	0.97478	0.97472	0.97466	0.97460	0.97454	0.97448	0.97442	0.97436	0.97430
80	0.97179	0.97173	0.97167	0.97160	0.97154	0.97148	0.97142	0.97135	0.97129	0.97123
85	0.96861	0.96855	0.96848	0.96842	0.96835	0.96829	0.96822	0.96816	0.96809	0.96803
90	0.96531	0.96524	0.96518	0.96511	0.96504	0.96497	0.96491	0.96484	0.96477	0.96470
95	0.96189	0.96182	0.96175	0.96168	0.96161	0.96154	0.96147	0.96140	0.96133	0.96126
"""

    # Medium pressure water density table (1,000 psi or 68.1 atm)
    medium_pressure_water_table = """
Temperature	Density (at 1000 psi or 68.1 atm)	Specific weight						
(°C)	(°F)	(g/cm3)	(kg/m3)	(sl/ft3)	(lbm/ft3)	(lbm/gal(US liq))	(lbf/ft3)	(kN/m3)
0.0	32	1.0031	1003.1	1.946	62.62	8.371	62.62	9.837
4.4	40	1.0031	1003.1	1.946	62.62	8.371	62.62	9.837
10.0	50	1.0031	1003.1	1.946	62.62	8.371	62.62	9.837
15.6	60	1.0024	1002.4	1.945	62.58	8.366	62.58	9.831
21.1	70	1.0012	1001.2	1.943	62.50	8.355	62.50	9.818
26.7	80	0.9999	999.9	1.940	62.42	8.344	62.42	9.805
32.2	90	0.9981	998.1	1.937	62.31	8.330	62.31	9.788
37.8	100	0.9962	996.2	1.933	62.19	8.314	62.19	9.769
43.3	110	0.9944	994.4	1.928	62.03	8.292	62.03	9.744
48.9	120	0.9912	991.2	1.923	61.88	8.272	61.88	9.721
54.4	130	0.9888	988.8	1.919	61.73	8.252	61.73	9.697
60.0	140	0.9864	986.4	1.914	61.58	8.232	61.58	9.673
65.6	150	0.9834	983.4	1.908	61.39	8.207	61.39	9.644
71.1	160	0.9803	980.3	1.902	61.20	8.181	61.20	9.614
76.7	170	0.9768	976.8	1.895	60.98	8.152	60.98	9.579
82.2	180	0.9731	973.1	1.888	60.75	8.121	60.75	9.543
87.8	190	0.9696	969.6	1.881	60.53	8.092	60.53	9.509
93.3	200	0.9661	966.1	1.874	60.31	8.062	60.31	9.474
121.1	250	0.9456	945.6	1.835	59.03	7.891	59.03	9.273
148.9	300	0.9217	921.7	1.788	57.54	7.692	57.54	9.039
176.7	350	0.8943	894.3	1.735	55.83	7.463	55.83	8.770
204.4	400	0.8636	863.6	1.676	53.91	7.207	53.91	8.469
260.0	500	0.7867	786.7	1.526	49.11	6.565	49.11	7.715
284.8	544.58	boiling point
"""

    # High pressure water density table (10,000 psi or 681 atm)
    high_pressure_water_table = """
Temperature	Density (at 10 000 psi or 681 atm)	Specific weight						
(°C)	(°F)	(g/cm3)	(kg/m3)	(sl/ft3)	(lbm/ft3)	(lbm/gal(US liq))	(lbf/ft3)	(kN/m3)
0.0	32	1.033	1033	2.004	64.5	8.62	64.5	10.13
4.4	40	1.032	1032	2.003	64.4	8.61	64.4	10.12
10.0	50	1.031	1031	2.000	64.4	8.60	64.4	10.11
15.6	60	1.029	1029	1.997	64.3	8.59	64.3	10.09
21.1	70	1.028	1028	1.994	64.1	8.58	64.1	10.08
26.7	80	1.026	1026	1.990	64.0	8.56	64.0	10.06
32.2	90	1.024	1024	1.986	63.9	8.54	63.9	10.04
37.8	100	1.021	1021	1.982	63.8	8.52	63.8	10.02
43.3	110	1.019	1019	1.977	63.6	8.51	63.6	9.99
48.9	120	1.017	1017	1.973	63.5	8.48	63.5	9.97
54.4	130	1.014	1014	1.968	63.3	8.46	63.3	9.94
60.0	140	1.011	1011	1.962	63.1	8.44	63.1	9.92
65.6	150	1.008	1008	1.957	63.0	8.42	63.0	9.89
71.1	160	1.005	1005	1.951	62.8	8.39	62.8	9.86
76.7	170	1.002	1002	1.945	62.6	8.37	62.6	9.83
82.2	180	0.999	999	1.939	62.4	8.34	62.4	9.80
87.8	190	0.996	996	1.932	62.2	8.31	62.2	9.77
93.3	200	0.992	992	1.926	62.0	8.28	62.0	9.73
121.1	250	0.974	974	1.890	60.8	8.13	60.8	9.55
148.9	300	0.953	953	1.849	59.5	7.95	59.5	9.35
176.7	350	0.930	930	1.805	58.1	7.76	58.1	9.12
204.4	400	0.905	905	1.756	56.5	7.55	56.5	8.88
260.0	500	0.847	847	1.643	52.9	7.07	52.9	8.31
315.6	600	0.774	774	1.501	48.3	6.46	48.3	7.59
"""
    
    # Parse the water density table data
    atm_temps_c, atm_densities = parse_atmospheric_water_table(atmospheric_water_table)
    med_temps_c, med_densities = parse_pressure_water_table(medium_pressure_water_table)
    high_temps_c, high_densities = parse_pressure_water_table(high_pressure_water_table)
    
    # Define pressure values in MPa
    atm_pressure = 0.101325  # Atmospheric pressure in MPa
    med_pressure = 6.9       # Medium pressure (1,000 psi or 68.1 atm) in MPa
    high_pressure = 69.0     # High pressure (10,000 psi or 681 atm) in MPa
    
    # Test API at each pressure level
    print("Testing water density API at atmospheric pressure...")
    atm_results = test_water_api(atm_temps_c, atm_densities, atm_pressure, 100.0)
    
    print("\nTesting water density API at medium pressure (6.9 MPa)...")
    med_results = test_water_api(med_temps_c, med_densities, med_pressure, 260.0)
    
    print("\nTesting water density API at high pressure (69 MPa)...")
    high_results = test_water_api(high_temps_c, high_densities, high_pressure, 315.6)
    
    # Create combined visualization only if we have valid results for all three pressures
    if atm_results and med_results and high_results:
        plt.figure(figsize=(15, 10))
        
        # Define color scheme for each pressure level
        colors = {
            'atmospheric': {'ref': 'blue', 'api': 'lightblue'},
            'medium': {'ref': 'green', 'api': 'lightgreen'},
            'high': {'ref': 'red', 'api': 'salmon'}
        }
        
        # Plot atmospheric pressure data
        plt.plot(atm_results['temps_c'], atm_results['ref_densities'], 
                 color=colors['atmospheric']['ref'], linewidth=2.5, 
                 label='Reference (Atmospheric)')
        plt.plot(atm_results['temps_c'], atm_results['api_densities'], 
                 color=colors['atmospheric']['api'], linewidth=2.0, linestyle='--', 
                 label='API (Atmospheric)')
        
        # Plot medium pressure data
        plt.plot(med_results['temps_c'], med_results['ref_densities'], 
                 color=colors['medium']['ref'], linewidth=2.5, 
                 label='Reference (1,000 psi)')
        plt.plot(med_results['temps_c'], med_results['api_densities'], 
                 color=colors['medium']['api'], linewidth=2.0, linestyle='--', 
                 label='API (1,000 psi)')
        
        # Plot high pressure data
        plt.plot(high_results['temps_c'], high_results['ref_densities'], 
                 color=colors['high']['ref'], linewidth=2.5, 
                 label='Reference (10,000 psi)')
        plt.plot(high_results['temps_c'], high_results['api_densities'], 
                 color=colors['high']['api'], linewidth=2.0, linestyle='--', 
                 label='API (10,000 psi)')
        
        # Add informative labels and title
        plt.xlabel('Temperature (°C)', fontsize=14)
        plt.ylabel('Density (kg/m³)', fontsize=14)
        plt.title('Water Density at Different Pressures: Calculated vs Reference Values', 
                  fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Improve axis readability
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        # Add minor gridlines for better readability
        plt.minorticks_on()
        plt.grid(which='minor', linestyle=':', alpha=0.2)
        
        # Add some padding around the plot
        plt.tight_layout(pad=2.0)
        
        # Save the figure with high resolution
        plt.savefig('combined_water_density_comparison.png', dpi=300)
        
        print("\nCombined visualization created and saved as 'combined_water_density_comparison.png'")
    else:
        print("\nUnable to create combined visualization - missing valid data for one or more pressure levels")
        
    # Print combined analysis
    print("\nOverall Analysis:")
    all_valid = True
    
    if atm_results and len(atm_results['diff']) > 0:
        atm_mean_diff = np.mean(atm_results['diff'])
        print(f"Average difference at atmospheric pressure: {atm_mean_diff:.2f} kg/m³")
    else:
        all_valid = False
        
    if med_results and len(med_results['diff']) > 0:
        med_mean_diff = np.mean(med_results['diff'])
        print(f"Average difference at medium pressure (6.9 MPa): {med_mean_diff:.2f} kg/m³")
    else:
        all_valid = False
        
    if high_results and len(high_results['diff']) > 0:
        high_mean_diff = np.mean(high_results['diff'])
        print(f"Average difference at high pressure (69 MPa): {high_mean_diff:.2f} kg/m³")
    else:
        all_valid = False
        
    if not all_valid:
        print("Unable to calculate average differences for all pressure levels - some data points missing")

if __name__ == "__main__":
    main()