# Calcium-Imaging-Analysis-DRD1-vs-DRD2
This project contains MATLAB code for analyzing and visualizing calcium imaging data recorded from DRD1 and DRD2 neurons and analyzed with Suite2P. The pipeline processes `.mat` files containing neuronal activity data, normalizes and clusters it, and generates informative heatmaps and plots.

Calcium_events_extraction.m # Main MATLAB script for data analysis
1. **Loads** all calcium traces from `.mat` files in DRD1 and DRD2 folders.
2. **Normalizes** each neuron's trace using baseline.
3. **Extract calcium transients**
   
Data_visualization.m # MATLAB script for visualization
1. **Sorts** neurons by early vs late activity change.
2. **Performs clustering** (k-means and PCA) to group neurons into functional types.
3. **Generates visualizations**:
    - Heatmaps 
    - Mean cluster traces
