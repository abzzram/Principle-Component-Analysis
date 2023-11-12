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
```
% which wells to input into the PCA
wells = [F.MCF10_EKAR35_EGF1 F.A1_EKAR35_EGF1 F.MCF7_EKAR35_EGF1 F.MCF10_EKAR35_AREG1 ...
    F.A1_EKAR35_AREG1 F.MCF7_EKAR35_AREG1 F.Imagingmedia]
% or, alternatively 
wells = [1:96]
```

Set input parameters for PCA in the v and v2 variables. 
```
v = {'samplingtime',6,'wells',wells,'chkpulse',0,'ploton',1}
% Pulse Analysis Definitions (see matalabs findpeaks() function for help)
v2 = {'narm',5,'min_rise',.06,'min_vy',0.06,'MAXW', 150}; 
```

For input definitions and other optional inputs, type:
```
help PCApl
```

Run the main PCA analysis using the PCApl function with the provided parameters.
```
[Z, PCA, PARAMS, Names]=PCApl(data,'ekar',datasheet,v,v2);
```

## Output
Three figures will pop up:
Figure 1: Percent varance explain by each PC
Figure 2: Weights of each component
Figure 3: 3D scatter plot of cells in PC space. There will be three scatter plots because they will each be colored by different categories