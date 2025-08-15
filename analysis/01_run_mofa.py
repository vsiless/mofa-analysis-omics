#!/usr/bin/env python3
"""
01_run_mofa.py
==============

Step 1: Run MOFA+ analysis on multi-omics data
This script loads the omics data and runs MOFA+ to extract latent factors.
"""

import pandas as pd
import numpy as np
from mofapy2.run.entry_point import entry_point
import os

def main():
    print("=" * 60)
    print("STEP 1: MOFA+ ANALYSIS")
    print("=" * 60)
    
    # Load CSVs
    print("Loading omics data...")
    lipidomics = pd.read_csv("omics/lipidomics.csv", index_col=0)
    metabolomics = pd.read_csv("omics/metabolomics.csv", index_col=0)
    proteomics = pd.read_csv("omics/proteomics.csv", index_col=0)
    transcriptomics = pd.read_csv("omics/transcriptomics.csv", index_col=0)

    # Convert to NumPy arrays
    data_arrays = [lipidomics.values, metabolomics.values, proteomics.values, transcriptomics.values]
    views_names = ["Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics"]
    likelihoods = ["gaussian"] * 4  # continuous omics

    # Sanity check: all arrays must have the same number of samples (rows)
    num_samples = [arr.shape[0] for arr in data_arrays]
    print(f"Number of samples per view: {dict(zip(views_names, num_samples))}")

    if len(set(num_samples)) != 1:
        raise ValueError("All views must have the same number of samples (rows). Check your CSVs!")

    # Initialize MOFA+
    print("Initializing MOFA+...")
    ent = entry_point()

    # Set data matrix - data should be a nested list: first dimension for views, second dimension for groups
    data_nested = [[arr] for arr in data_arrays]  # Each view is a list containing one group

    ent.set_data_matrix(
        data=data_nested,        # Nested list: views x groups
        views_names=views_names,   # flat list, one per view
        likelihoods=likelihoods    # flat list
    )

    # Set options
    print("Setting MOFA+ options...")
    ent.set_data_options(scale_views=True)
    ent.set_model_options(factors=10)
    ent.set_train_options(iter=1000, convergence_mode="fast")

    # Run MOFA+
    print("Running MOFA+ training...")
    ent.build()
    ent.run()
    
    # Save model
    print("Saving MOFA+ model...")
    os.makedirs("results",exist_ok=True)
    ent.save("results/mofa_model.hdf5")

    print("✅ MOFA+ analysis complete! Model saved as mofa_model.hdf5")
    print("=" * 60)

if __name__ == "__main__":
    main()
