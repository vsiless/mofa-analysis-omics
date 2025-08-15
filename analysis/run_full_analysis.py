#!/usr/bin/env python3
"""
run_full_analysis.py
===================

Master script to run the complete MOFA+ analysis pipeline
This script runs all three analysis steps in order:
1. Run MOFA+ analysis
2. Generate plots
3. Perform clustering analysis
"""

import subprocess
import sys
import os
from datetime import datetime

def run_script(script_name, description):
    """Run a Python script and handle errors"""
    print(f"\n{'='*60}")
    print(f"RUNNING: {script_name}")
    print(f"DESCRIPTION: {description}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            print("✅ SUCCESS!")
            print(result.stdout)
        else:
            print("❌ ERROR!")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ EXCEPTION: {e}")
        return False
    
    return True

def main():
    print("🚀 MOFA+ COMPLETE ANALYSIS PIPELINE")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Working directory: {os.getcwd()}")
    
    # Check if we're in the analysis directory
    if not os.path.exists("01_run_mofa.py"):
        print("❌ ERROR: Please run this script from the 'analysis' directory")
        print("   cd analysis")
        print("   python run_full_analysis.py")
        return
    
    # Step 1: Run MOFA+ analysis
    if not run_script("01_run_mofa.py", "Run MOFA+ analysis on multi-omics data"):
        print("❌ Step 1 failed. Stopping analysis.")
        return
    
    # Step 2: Generate plots
    if not run_script("02_plot_results.py", "Generate MOFA+ result plots"):
        print("❌ Step 2 failed. Stopping analysis.")
        return
    
    # Step 3: Clustering analysis
    if not run_script("03_clustering_analysis.py", "Perform clustering analysis"):
        print("❌ Step 3 failed. Stopping analysis.")
        return
    
    print(f"\n{'='*60}")
    print("🎉 COMPLETE ANALYSIS PIPELINE FINISHED!")
    print(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")
    
    print("\n📁 OUTPUT FILES:")
    print("   📊 MOFA+ Model: ../mofa_model.hdf5")
    print("   📈 Figures: ../figs/ (10 PNG files)")
    print("   📋 Clustering Results: ../mofa_clustering_results.csv")
    
    print("\n📊 FIGURES GENERATED:")
    print("   01_variance_explained_heatmap.png")
    print("   02_total_variance_explained.png")
    print("   03_factor_distributions.png")
    print("   04_factor_correlations.png")
    print("   05_feature_weights_lipidomics.png")
    print("   05_feature_weights_metabolomics.png")
    print("   05_feature_weights_proteomics.png")
    print("   06_elbo_convergence.png")
    print("   07_clustering_dendrogram.png")
    print("   08_pca_with_clusters.png")
    print("   09_factor_heatmap_by_cluster.png")
    print("   10_cluster_characterization.png")
    
    print("\n✅ Analysis complete! Check the figs/ folder for all visualizations.")

if __name__ == "__main__":
    main()
