# CellPainting

## Analysis of Cell Painting morphological data

### Requirements  
Aside from standard packages (numpy, pandas, sklearn, matplotlib):  
anndata (tested using 0.10.7)  
scanpy (testing using 1.10.1)

### Introductions  

This project was intended to understand how cell morphology changes across drug conditions using Cell Painting (CP). 
The analysis is performed at the well-level, in which the CP features are averaged across all cells within each well. 
  
The data for this analysis comes in three parts:
  - an excel file that maps each drug code to the plate well row and column
  - an excel file that maps the drug code to the drug name and concentration
  - a list of directories, each including a .txt file with the well-averaged CP features, and (optional) the well-averaged phagocytosis data.

These three data sources are found in the data directory, which was the result of a prelinary drug screen testing the effect of ~33 drugs on HMC3 phagocytosis.

In analysis.py, the class CreateData will create an h5ad (AnnData) file, in which .X matrix contains the normalized CP data, obs
contains the well location, drug name and concentration, and (optional) phagocytosis data, and var contains the CP feature names.
The class Analysis allows one to plot the various dose-response curves, perform PLS regression, etc. This file should be customized based on the user's desired goals.




