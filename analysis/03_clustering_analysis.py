#!/usr/bin/env python3
"""
03_clustering_analysis.py
=========================

Step 3: Perform clustering analysis on MOFA+ factors
This script performs hierarchical clustering on the factor values and generates
clustering visualizations. All figures are saved to the figs/ folder.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
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
    
    return results

def main():
    print("=" * 60)
    print("STEP 3: CLUSTERING ANALYSIS")
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
    view_names = ["Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics"]
    
    print(f"Data shapes: Factors {factors.shape}")
    
    # Prepare data for clustering
    print("Preparing data for clustering...")
    scaler = StandardScaler()
    factors_scaled = scaler.fit_transform(factors)
    print(f"Standardized factors - Mean: {factors_scaled.mean():.3f}, Std: {factors_scaled.std():.3f}")
    
    # Perform hierarchical clustering
    print("Performing hierarchical clustering...")
    linkage_matrix = linkage(factors_scaled, method='ward')
    
    # Plot dendrogram
    print("1. Generating clustering dendrogram...")
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, 
              labels=[f'Sample_{i+1}' for i in range(len(factors_scaled))],
              orientation='top',
              distance_sort='descending',
              show_leaf_counts=True)
    
    plt.title('Hierarchical Clustering Dendrogram', fontsize=16)
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig("figs/07_clustering_dendrogram.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Determine optimal number of clusters
    optimal_k = 3  # Based on typical biological patterns
    print(f"Using {optimal_k} clusters for analysis...")
    
    # Cut dendrogram to get clusters
    hierarchical_clusters = fcluster(linkage_matrix, optimal_k, criterion='maxclust')
    hierarchical_clusters = hierarchical_clusters - 1  # Convert to 0-based indexing
    
    print(f"Cluster sizes: {np.bincount(hierarchical_clusters)}")
    
    # 2. PCA plot of factors with clusters
    print("2. Generating PCA plot with clusters...")
    pca = PCA(n_components=2)
    factors_pca = pca.fit_transform(factors_scaled)
    
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(factors_pca[:, 0], factors_pca[:, 1], 
                         c=hierarchical_clusters, cmap='viridis', s=100, alpha=0.7)
    plt.colorbar(scatter, label='Cluster')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.title('PCA of Factor Values with Clusters', fontsize=16)
    plt.grid(True, alpha=0.3)
    
    # Add sample labels
    for i, (x, y) in enumerate(factors_pca):
        plt.annotate(f'S{i+1}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.tight_layout()
    plt.savefig("figs/08_pca_with_clusters.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Heatmap of factor values by cluster
    print("3. Generating factor values heatmap by cluster...")
    plt.figure(figsize=(12, 8))
    cluster_order = np.argsort(hierarchical_clusters)
    factor_data_ordered = factors_scaled[cluster_order]
    cluster_labels_ordered = hierarchical_clusters[cluster_order]
    
    sns.heatmap(factor_data_ordered.T, 
                cmap='RdBu_r', 
                center=0,
                xticklabels=[f'S{i+1}' for i in cluster_order],
                yticklabels=[f'Factor {i+1}' for i in range(factors.shape[1])],
                cbar_kws={'label': 'Standardized Factor Value'})
    
    # Add cluster boundaries
    cluster_boundaries = np.where(np.diff(cluster_labels_ordered) != 0)[0] + 0.5
    for boundary in cluster_boundaries:
        plt.axvline(x=boundary, color='black', linewidth=2)
    
    plt.title('Factor Values Heatmap by Cluster', fontsize=16)
    plt.xlabel('Samples (ordered by cluster)', fontsize=14)
    plt.ylabel('Factors', fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("figs/09_factor_heatmap_by_cluster.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Cluster characterization
    print("4. Generating cluster characterization plots...")
    cluster_means = []
    cluster_stds = []
    
    for cluster_id in range(optimal_k):
        cluster_mask = hierarchical_clusters == cluster_id
        cluster_data = factors_scaled[cluster_mask]
        
        cluster_means.append(np.mean(cluster_data, axis=0))
        cluster_stds.append(np.std(cluster_data, axis=0))
    
    cluster_means = np.array(cluster_means)
    cluster_stds = np.array(cluster_stds)
    
    # Plot cluster profiles
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Mean factor values per cluster
    im1 = ax1.imshow(cluster_means.T, cmap='RdBu_r', aspect='auto')
    ax1.set_xlabel('Cluster')
    ax1.set_ylabel('Factor')
    ax1.set_title('Mean Factor Values per Cluster')
    ax1.set_xticks(range(optimal_k))
    ax1.set_xticklabels([f'Cluster {i+1}' for i in range(optimal_k)])
    ax1.set_yticks(range(factors.shape[1]))
    ax1.set_yticklabels([f'Factor {i+1}' for i in range(factors.shape[1])])
    plt.colorbar(im1, ax=ax1, label='Standardized Factor Value')
    
    # Standard deviation per cluster
    im2 = ax2.imshow(cluster_stds.T, cmap='viridis', aspect='auto')
    ax2.set_xlabel('Cluster')
    ax2.set_ylabel('Factor')
    ax2.set_title('Factor Variability per Cluster (Std Dev)')
    ax2.set_xticks(range(optimal_k))
    ax2.set_xticklabels([f'Cluster {i+1}' for i in range(optimal_k)])
    ax2.set_yticks(range(factors.shape[1]))
    ax2.set_yticklabels([f'Factor {i+1}' for i in range(factors.shape[1])])
    plt.colorbar(im2, ax=ax2, label='Standard Deviation')
    
    plt.tight_layout()
    plt.savefig("figs/10_cluster_characterization.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create summary dataframe
    print("5. Creating clustering results summary...")
    summary_data = []
    for i in range(len(factors)):
        row = {
            'Sample': f'Sample_{i+1}',
            'Cluster': hierarchical_clusters[i] + 1
        }
        
        # Add factor values
        for j in range(factors.shape[1]):
            row[f'Factor_{j+1}_Raw'] = factors[i, j]
            row[f'Factor_{j+1}_Scaled'] = factors_scaled[i, j]
        
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save results
    summary_df.to_csv("results/mofa_clustering_results.csv", index=False)
    print("   Results saved to: results/mofa_clustering_results.csv")
    
    # Print cluster characteristics
    print(f"\n📊 CLUSTER CHARACTERISTICS:")
    for cluster_id in range(optimal_k):
        cluster_mask = hierarchical_clusters == cluster_id
        cluster_size = np.sum(cluster_mask)
        cluster_samples = [f'S{i+1}' for i in np.where(cluster_mask)[0]]
        
        print(f"\nCluster {cluster_id + 1} (n={cluster_size}):")
        print(f"  Samples: {', '.join(cluster_samples)}")
        
        # Find most distinctive factors
        cluster_mean = cluster_means[cluster_id]
        abs_means = np.abs(cluster_mean)
        top_factors = np.argsort(abs_means)[-3:][::-1]
        
        print(f"  Most distinctive factors:")
        for factor_idx in top_factors:
            value = cluster_mean[factor_idx]
            print(f"    Factor {factor_idx + 1}: {value:+.2f}")
    
    # Final summary statistics
    print(f"\n✅ CLUSTERING ANALYSIS COMPLETE!")
    print(f"📁 Figures saved to: ../figs/")
    print(f"📊 Final Summary:")
    print(f"  Total samples: {len(factors)}")
    print(f"  Number of clusters: {optimal_k}")
    print(f"  Clustering method: Hierarchical (Ward)")
    
    print(f"\nCluster sizes:")
    for i in range(optimal_k):
        size = np.sum(hierarchical_clusters == i)
        print(f"  Cluster {i+1}: {size} samples ({size/len(factors):.1%})")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
