# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:04:19 2025

@author: dtf8829
"""

import os
import numpy as np
import pandas as pd
from fooof import FOOOFGroup

# --------- config ---------
folder_path = r'F:/testPSDs2'
complete_dir = os.path.join(folder_path, 'complete')
os.makedirs(complete_dir, exist_ok=True)

f_range = (1.0, 30.0)           # Hz
batch_size = 10_000             # process this many files per batch
n_jobs = max(1, os.cpu_count()-1)
# --------------------------

# Build list of files to process
all_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
done = set(os.listdir(complete_dir))
files = [f for f in all_files
         if f not in done and 'log' not in f and 'rel' not in f]

print(f"Found {len(files)} files to process.")

def load_two_cols(fp):
    # Fast, minimal read: first two columns only, float32
    return pd.read_csv(fp, usecols=[0,1], dtype=np.float32, engine='c', memory_map=True)

for start in range(0, len(files), batch_size):
    batch_files = files[start:start+batch_size]
    print(f"\nBatch {start//batch_size + 1}: {len(batch_files)} files")

    # --- Build frequency mask from the first file in the batch ---
    first_df = load_two_cols(os.path.join(folder_path, batch_files[0]))
    frex = first_df.iloc[:,1].to_numpy()
    mask = (frex >= f_range[0]) & (frex <= f_range[1])

    # Preallocate power matrix (n_spectra x n_freqs_in_range)
    nF = int(mask.sum())
    power_mat = np.empty((len(batch_files), nF), dtype=np.float32)

    # Load batch
    bad_idx = []  # keep track of unusable spectra
    for i, fn in enumerate(batch_files):
        if i % 1000 == 0: 
            print(f"............. {i}")
        try:
            df = load_two_cols(os.path.join(folder_path, fn))
            # If frequency axis doesn’t match exactly, you could realign here;
            # this version assumes identical freqs within a batch.
            pwr = df.iloc[:,0].to_numpy()
            power_mat[i, :] = pwr[mask]
        except Exception as e:
            bad_idx.append(i)
            power_mat[i, :] = np.nan
    print(".........spectra loaded")
    # Drop bad rows before fitting
    if bad_idx:
        keep = np.setdiff1d(np.arange(len(batch_files)), np.array(bad_idx))
        batch_files_fit = [batch_files[i] for i in keep]
        power_fit = power_mat[keep, :]
    else:
        batch_files_fit = batch_files
        power_fit = power_mat

    # Guard: skip empty batch
    if power_fit.size == 0:
        print("All spectra in this batch failed to load; skipping.")
        continue

    frex_sel = frex[mask]

    # --- Fit group in parallel ---
    fg = FOOOFGroup(peak_width_limits=[1.0, 8.0],
                    max_n_peaks=6,
                    min_peak_height=0.1,
                    peak_threshold=2.0,
                    aperiodic_mode='fixed',
                    verbose=False)
    fg.fit(frex_sel, power_fit, f_range, n_jobs=n_jobs)

    # Extract aperiodic params (offset, exponent) for each spectrum
    ap = fg.get_params('aperiodic_params')  # shape: (n_spectra_fit, 2)
    print(".........spectra fit")

    # --- Save individual outputs per file ---
    for i, fn in enumerate(batch_files_fit):
        try:
            src = os.path.join(folder_path, fn)
            dst = os.path.join(complete_dir, fn)

            # Re-read original CSV to preserve all columns,
            # then append the two FOOOF params.
            data_full = pd.read_csv(src, engine='c')
            data_full['offset'] = ap[i, 0]
            data_full['exponent'] = ap[i, 1]
            data_full.to_csv(dst, index=False)
        except Exception as e:
            print(f"Save failed for {fn}: {e}")

    # Optionally, mark the failures too (so they don’t retry next run)
    for i in bad_idx:
        fn = batch_files[i]
        # Create a tiny flag file to mark as attempted
        open(os.path.join(complete_dir, fn + '.failed'), 'a').close()

print("\nDone.")
