# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:25:38 2024

@author: rens_
"""
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import pandas as pd
from utils import print_underlined

smps_bins = [7.0, 11.5, 15.4, 20.5, 27.4, 36.5, 48.7, 64.9, 86.6, 115.5, 154.0, 205.4, 273.8, 365.2]
smps_counts = np.array([3.21707011e+08, 6.71035657e+08, 6.84832623e+08, 9.28110838e+08,
       9.99454158e+08, 9.12587603e+08, 7.68542700e+08, 6.23995149e+08,
       4.22267614e+08, 1.80658013e+08, 2.97094689e+07, 2.37900754e+06,
       8.96179783e+05])

salsa_bins = [3., 7.66309432, 19.57433821, 50., 96.71682101, 187.08286934, 361.88120777]

def translate_particles(original_bins, original_counts, new_bins):
    # Calculate the area under the curve for the original counts
    new_bin_sizes = np.diff(new_bins)
    area_original = trapz(original_counts, x=original_bins[:-1])
    # area_original = sum(np.diff(original_bins) * original_counts)

    # Initialize array to store translated counts for new bins
    translated_counts = np.zeros(len(new_bins) - 1)

    # Translate particle counts to new bins
    for i in range(len(new_bins) - 1):
        # Calculate the overlap between the original and new bins
        overlap = np.maximum(0, np.minimum(original_bins[1:], new_bins[i+1]) - np.maximum(original_bins[:-1], new_bins[i]))
        overlap_ratio = overlap / new_bin_sizes[i]

        # Allocate particles from each original bin to the new bin
        translated_counts[i] = np.sum(overlap_ratio * original_counts)

    # Calculate the actual area under the curve for the translated counts
    area_translated = trapz(translated_counts, x=new_bins[:-1])
    # area_translated = sum(np.diff(new_bins) * translated_counts)

    # Scale the translated counts to match the total count
    translated_counts /= (area_original / area_translated)
    # test = pd.DataFrame(data = [translated_counts], columns = ["1a1", "1a2", "1a3", "2a1", "2a2", "2a3"])

    return translated_counts

smps_translated = translate_particles(smps_bins, smps_counts, salsa_bins)

plt.plot(salsa_bins[1:], smps_translated, label = "SALSA2.0")
plt.plot(smps_bins[1:], smps_counts, label = "SMPS")
plt.legend()
plt.title("Translation from SMPS to SALSA bins")
plt.xlabel("Bin boundaries (nm)")
plt.ylabel("# particles m$^{-3}$")

for i in range(len(salsa_bins) - 1):
    plt.fill_between(salsa_bins[i+1:i+3], smps_translated[i:i+2], 0, color = "blue", alpha = 0.3)
    
plt.show()

print_underlined("\nBy (rectanguar) integration")
print(f"Total in original: {sum(np.diff(smps_bins) * smps_counts):.2e}")
print(f"Total in translated: {sum(np.diff(salsa_bins) * smps_translated):.2e}")
print(f"Ratio: {sum(np.diff(salsa_bins) * smps_translated)/sum(np.diff(smps_bins) * smps_counts):.2f}")

print_underlined("\nBy (trapezoidal) integration")
print(f"Total in original: {trapz(smps_counts, smps_bins[1:]):.2e}")
print(f"Total in translated: {trapz(smps_translated, salsa_bins[1:]):.2e}")
print(f"Ratio: {trapz(smps_translated, salsa_bins[1:])/trapz(smps_counts, smps_bins[1:]):.2f}")
