#!/usr/bin/env python3
"""
02_plot_results.py
==================

Step 2: Generate MOFA+ result plots
This script loads the MOFA+ model and generates all visualization plots.
All figures are saved to the figs/ folder.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import os
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

def load_mofa_results(model_path="results/mofa_model.hdf5"):
    """Load results from the saved MOFA+ model"""
    results = {}
    
    with h5py.File(model_path, 'r') as f:
        # Load factor values (Z) - organized by groups
        if 'expectations' in f and 'Z' in f['expectations']:
            group_name = list(f['expectations']['Z'].keys())[0]
            results['factors'] = f['expectations']['Z'][group_name][:].T  # samples x factors
        
        # Load weights (W) - organized by views
        if 'expectations' in f and 'W' in f['expectations']:
            weights = []
            view_names = ["Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics"]
            for view in view_names:
                if view in f['expectations']['W']:
                    weights.append(f['expectations']['W'][view][:])
            results['weights'] = weights
        
        # Load variance explained
        if 'variance_explained' in f:
            group_name = list(f['variance_explained']['r2_per_factor'].keys())[0]
            results['variance_explained'] = f['variance_explained']['r2_per_factor'][group_name][:]
        
        # Load training stats
        if 'training_stats' in f:
            results['training_stats'] = {}
            for key in f['training_stats'].keys():
                if hasattr(f['training_stats'][key], '__getitem__'):
                    results['training_stats'][key] = f['training_stats'][key][:]
    
    return results

def main():
    print("=" * 60)
    print("STEP 2: GENERATE MOFA+ PLOTS")
    print("=" * 60)
    
    # Create figs directory if it doesn't exist
    os.makedirs("figs", exist_ok=True)
    
    # Load the results
    print("Loading MOFA+ results...")
    results = load_mofa_results()
    print("Results loaded successfully!")
    
    # Extract key data
    factors = results['factors']
    var_explained = results['variance_explained']
    weights = results['weights']
    view_names = ["Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics"]
    
    print(f"Data shapes: Factors {factors.shape}, Variance explained {var_explained.shape}")
    
    # 1. Variance explained heatmap
    print("1. Generating variance explained heatmap...")
    plt.figure(figsize=(14, 8))
    sns.heatmap(var_explained, 
                annot=True, 
                fmt='.1f', 
                cmap='YlOrRd',
                xticklabels=[f'Factor {i+1}' for i in range(var_explained.shape[1])],
                yticklabels=view_names,
                cbar_kws={'label': 'Variance Explained (%)'})
    plt.title('Variance Explained by Each Factor Across Views', fontsize=16, pad=20)
    plt.xlabel('Factors', fontsize=14)
    plt.ylabel('Views', fontsize=14)
    plt.tight_layout()
    plt.savefig("figs/01_variance_explained_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Total variance explained per view
    print("2. Generating total variance explained plot...")
    total_var_per_view = np.sum(var_explained, axis=1)
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(view_names, total_var_per_view, color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
    
    # Add value labels on bars
    for bar, value in zip(bars, total_var_per_view):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    plt.title('Total Variance Explained by All Factors per View', fontsize=16, pad=20)
    plt.ylabel('Variance Explained (%)', fontsize=14)
    plt.xlabel('Views', fontsize=14)
    plt.ylim(0, max(total_var_per_view) * 1.1)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig("figs/02_total_variance_explained.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Factor distributions
    print("3. Generating factor distributions...")
    n_factors = factors.shape[1]
    n_cols = 3
    n_rows = (n_factors + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_factors > 1 else [axes]
    
    for i in range(n_factors):
        ax = axes[i]
        factor_values = factors[:, i]
        
        # Create histogram
        ax.hist(factor_values, bins=15, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_title(f'Factor {i+1} Distribution', fontsize=12)
        ax.set_xlabel('Factor Value')
        ax.set_ylabel('Frequency')
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for i in range(n_factors, len(axes)):
        axes[i].set_visible(False)
    
    plt.suptitle('Distribution of Factor Values Across Samples', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig("figs/03_factor_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Factor correlations
    print("4. Generating factor correlation matrix...")
    corr_matrix = np.corrcoef(factors.T)
    
    plt.figure(figsize=(10, 8))
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, 
                mask=mask,
                annot=True, 
                fmt='.2f', 
                cmap='RdBu_r',
                center=0,
                square=True,
                xticklabels=[f'Factor {i+1}' for i in range(corr_matrix.shape[0])],
                yticklabels=[f'Factor {i+1}' for i in range(corr_matrix.shape[0])],
                cbar_kws={'label': 'Correlation Coefficient'})
    plt.title('Factor Correlation Matrix', fontsize=16, pad=20)
    plt.xlabel('Factors', fontsize=14)
    plt.ylabel('Factors', fontsize=14)
    plt.tight_layout()
    plt.savefig("figs/04_factor_correlations.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Feature weights for each view (excluding transcriptomics due to size)
    print("5. Generating feature weights plots...")
    for view_idx, view_name in enumerate(["Lipidomics", "Metabolomics", "Proteomics"]):
        print(f"   Plotting {view_name} feature weights...")
        view_weights = weights[view_idx]
        
        # Create subplots for each factor
        n_factors = view_weights.shape[1]
        n_cols = 2
        n_rows = (n_factors + n_cols - 1) // n_cols
        
        # Adjust figure size to prevent "too large" error
        fig_height = min(6 * n_rows, 30)  # Cap height at 30 inches
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, fig_height))
        axes = axes.flatten() if n_factors > 1 else [axes]
        
        for i in range(n_factors):
            ax = axes[i]
            factor_weights = view_weights[:, i]
            
            # Get top features (absolute values)
            abs_weights = np.abs(factor_weights)
            top_n = min(20, len(factor_weights))
            top_indices = np.argsort(abs_weights)[-top_n:]
            
            # Plot top features
            y_pos = np.arange(len(top_indices))
            ax.barh(y_pos, factor_weights[top_indices], color='lightcoral')
            ax.set_yticks(y_pos)
            ax.set_yticklabels([f'Feature {idx+1}' for idx in top_indices])
            ax.set_title(f'{view_name} - Factor {i+1} Top {top_n} Features', fontsize=12)
            ax.set_xlabel('Weight')
            ax.grid(True, alpha=0.3)
        
        # Hide empty subplots
        for i in range(n_factors, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Top Feature Weights for {view_name}', fontsize=16, y=0.98)
        plt.tight_layout()
        plt.savefig(f"figs/05_feature_weights_{view_name.lower()}.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. ELBO convergence (if available)
    print("6. Generating ELBO convergence plot...")
    if 'training_stats' in results and 'elbo' in results['training_stats']:
        elbo_values = results['training_stats']['elbo']
        
        plt.figure(figsize=(10, 6))
        plt.plot(elbo_values, 'b-o', linewidth=2, markersize=6)
        plt.title('ELBO Convergence During Training', fontsize=16, pad=20)
        plt.xlabel('Iteration', fontsize=12)
        plt.ylabel('ELBO', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig("figs/06_elbo_convergence.png", dpi=300, bbox_inches='tight')
        plt.close()
    else:
        print("   ELBO convergence data not available")
    
    # Print summary
    print(f"\n✅ All plots generated successfully!")
    print(f"📁 Figures saved to: ../figs/")
    print(f"📊 Total variance explained per view:")
    for i, view in enumerate(view_names):
        print(f"   {view}: {total_var_per_view[i]:.1f}%")
    print("=" * 60)

if __name__ == "__main__":
    main()
