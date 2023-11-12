# Principle Component Analysis
 These functions are created for standard analysis of single-cell biosensor data.

These functions will featurize the single-cell time series traces into several parameters of biosensor activity.

## Getting Started
To run the PCA analysis, follow these steps:

##  Prerequisites
Make sure you have MATLAB installed on your system, prefferably 2020 or newer. 

## Usage
Open MATLAB and load your dataset.
``` 
data = load('example_data.mat','dcat')

```

```
% define location of the experiment datasheet 
datasheet = 'L:\Databases\ImagingExperimentSheets\Ram experiments\2020-06-0210A_184_7_EGF_EREG.xlsx'
```
Optional (Run the function datasheetIndexAR_Tx1 to extract cell type and treatment labels.) If you decide to not use this, manually input the well # into PCApl

```
% extract cell type and treatment labels for wells in the data
[F,L] = datasheetIndexAR_Tx1(datasheet);


```
Define the wells for PCA by specifying the wells variable.
Set input parameters for PCA in the v and v2 variables.
Run the main PCA analysis using the PCApl function with the provided parameters.