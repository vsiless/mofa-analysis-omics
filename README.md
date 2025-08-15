# MOFA+ Multi-Omics Analysis Pipeline

This folder contains a clean, organized pipeline for MOFA+ multi-omics analysis with clustering.

## 📁 Folder Structure

```
mofa/
├── omics/                    # Input data
│   ├── lipidomics.csv
│   ├── metabolomics.csv
│   ├── proteomics.csv
│   └── transcriptomics.csv
├── analysis/                 # Analysis scripts (this folder)
│   ├── 01_run_mofa.py       # Step 1: Run MOFA+ analysis
│   ├── 02_plot_results.py   # Step 2: Generate plots
│   ├── 03_clustering_analysis.py  # Step 3: Clustering analysis
│   ├── run_full_analysis.py # Master script to run everything
│   └── README.md            # This file
├── figs/                     # Output figures (created automatically)
└── mofa_model.hdf5          # MOFA+ model (created by Step 1)
```

## 🚀 Quick Start

### Option 1: Run Complete Pipeline
```bash
cd analysis
python run_full_analysis.py
```

### Option 2: Run Steps Individually
```bash
cd analysis

# Step 1: Run MOFA+ analysis
python 01_run_mofa.py

# Step 2: Generate plots
python 02_plot_results.py

# Step 3: Perform clustering
python 03_clustering_analysis.py
```

## 📊 Analysis Steps

### Step 1: `01_run_mofa.py`
- Loads multi-omics data from CSV files
- Runs MOFA+ analysis with 10 factors
- Saves model as `mofa_model.hdf5`

### Step 2: `02_plot_results.py`
- Loads MOFA+ model results
- Generates 6 visualization plots:
  - Variance explained heatmap
  - Total variance explained per view
  - Factor distributions
  - Factor correlations
  - Feature weights (for each view)
  - ELBO convergence (if available)
- All figures saved to `figs/` folder

### Step 3: `03_clustering_analysis.py`
- Performs hierarchical clustering on factor values
- Generates 4 clustering visualizations:
  - Clustering dendrogram
  - PCA plot with clusters
  - Factor heatmap by cluster
  - Cluster characterization
- Saves clustering results as CSV
- All figures saved to `figs/` folder

## 📈 Output Files

### Figures (in `figs/` folder)
1. `01_variance_explained_heatmap.png` - Variance explained by each factor
2. `02_total_variance_explained.png` - Total variance per view
3. `03_factor_distributions.png` - Distribution of factor values
4. `04_factor_correlations.png` - Factor correlation matrix
5. `05_feature_weights_*.png` - Top feature weights per view
6. `06_elbo_convergence.png` - Training convergence
7. `07_clustering_dendrogram.png` - Hierarchical clustering tree
8. `08_pca_with_clusters.png` - PCA plot with cluster colors
9. `09_factor_heatmap_by_cluster.png` - Factor values by cluster
10. `10_cluster_characterization.png` - Cluster profiles

### Data Files
- `mofa_model.hdf5` - Trained MOFA+ model
- `mofa_clustering_results.csv` - Sample clustering assignments

## 🔧 Requirements

Make sure you have the following packages installed:
```bash
pip install mofapy2 pandas numpy matplotlib seaborn h5py scikit-learn scipy
```

## 📋 Expected Results

Based on your data, you should expect:
- **24 samples** across 4 omics views
- **10 latent factors** extracted by MOFA+
- **3 clusters** identified by hierarchical clustering
- **High variance explained** in metabolomics and proteomics
- **Lower variance explained** in transcriptomics

## 🎯 Key Insights

The analysis will reveal:
1. **Variance explained** by each factor across omics views
2. **Factor correlations** and independence
3. **Sample clusters** based on multi-omics profiles
4. **Most important features** contributing to each factor
5. **Biological patterns** captured by the latent factors

## 🐛 Troubleshooting

- **Missing files**: Make sure all CSV files are in the `omics/` folder
- **Import errors**: Install required packages with pip
- **Memory issues**: Reduce number of factors in `01_run_mofa.py`
- **Plot errors**: Check matplotlib backend settings

## 📝 Customization

You can modify the analysis by editing:
- **Number of factors**: Change `factors=10` in `01_run_mofa.py`
- **Number of clusters**: Change `optimal_k = 3` in `03_clustering_analysis.py`
- **Plot styles**: Modify matplotlib settings in each script
- **Output formats**: Change file extensions and DPI settings
