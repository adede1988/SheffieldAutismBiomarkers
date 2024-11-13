# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:10:10 2024

@author: Adam Dede
"""

import os
import pandas as pd
import numpy as np
from fooof import FOOOF
import matplotlib.pyplot as plt
from scipy.stats import pearsonr 

# Define the directory containing your CSV files
folder_path = 'H:/testPSDs'  # Update this to your folder path

# Create an empty list to store exponents and offsets
results = []

# Define the frequency range to fit
f_range = [1, 30]  # Adjust as needed for your use case

fileList = os.listdir(folder_path)

# Loop through all CSV files in the folder
for ii in range(1715):
    if fileList[ii].endswith('.csv'):
        # Construct the full path to the CSV file
        file_path = os.path.join(folder_path, fileList[ii])

        # Read the CSV file using pandas
        data = pd.read_csv(file_path)

        # Assuming the first column is power and the second column is frequency
        power_values = data.iloc[:, 0].values  # Extract power values from the first column
        frequency_values = data.iloc[:, 1].values  # Extract frequency values from the second column

        # Convert frequency and power values to NumPy arrays
        frex = np.array(frequency_values)
        psd = np.array(power_values)

        # Initialize the FOOOF model
        fooof_obj = FOOOF(peak_width_limits=[1.0, 8.0], max_n_peaks=6, min_peak_height=0.1,
                                     peak_threshold=2.0, aperiodic_mode='fixed')

        # Fit the FOOOF model
        fooof_obj.fit(frex, psd, f_range)

        # Get FOOOF results
        exponent = fooof_obj.aperiodic_params_[1]  # FOOOF exponent (slope)
        offset = fooof_obj.aperiodic_params_[0]    # FOOOF offset
        
        # Append the exponent and offset to the results list
        results.append([exponent, offset, data.iloc[1,2], data.iloc[1,3]])
        
        
# Convert the results list to a NumPy array for easier manipulation
results_array = np.array(results)

# Extract exponents and offsets for plotting
exponents = results_array[:, 0]
exponentsOG = results_array[:, 2]

# Create a scatter plot of exponents vs. offsets
plt.scatter(exponents, exponentsOG)
plt.title('slope method comparison')
plt.xlabel('Fooof calculated value')
plt.ylabel('Voytek 2015 calcualted value')
plt.grid(True)
plt.show()

# Perform Pearson correlation test
correlation_coefficient, p_value = pearsonr(exponents, exponentsOG)

# Print the correlation results
print(f"Pearson Correlation Coefficient: {correlation_coefficient}")
print(f"P-value: {p_value}")


# Extract exponents and offsets for plotting
exponents = results_array[:, 1]
exponentsOG = results_array[:, 3]

# Create a scatter plot of exponents vs. offsets
plt.scatter(exponents, exponentsOG)
plt.title('offset method comparison')
plt.xlabel('Fooof calculated value')
plt.ylabel('Voytek 2015 calcualted value')
plt.grid(True)
plt.show()

# Perform Pearson correlation test
correlation_coefficient, p_value = pearsonr(exponents, exponentsOG)

# Print the correlation results
print(f"Pearson Correlation Coefficient: {correlation_coefficient}")
print(f"P-value: {p_value}")

