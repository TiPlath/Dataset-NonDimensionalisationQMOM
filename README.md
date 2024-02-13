*Version: 1.0, Date: 16.02.2023*
 
### __Data to reproduce the paper: "Non-dimensionalization of quadrature method of moments for wet granulation"__

__Authors: Plath, Timo¹*; Luding, Stefan¹; Weinhart, Thomas¹__

¹ Thermal and Fluid Engineering (ET) | University of Twente | P.O. Box 217, 7500 AE Enschede, The Netherlands |  
\* Corresponding author: Plath, Timo; t.plath@utwente.nl

***General Introduction:***
The software consists of a fully working MATLAB (Version R2019b) script for a non-dimensional and dimensional quadrature method of moments
and all the additional methods that are needed to run the quadrature method of moments (reading in data, adaptive Wheeler algorithm, etc.).
Additionally a Python script (Version 3.7) for the Buckingham theorem using BuckinghamPy is available.
Furthermore this set of software allows to reproduce all figures and verifies integrity of the title paper. It makes the software reusable and accessible according to the FAIR principles.  
This research was funded by NWO VIDI grant 16604, *“Virtual Prototyping of Particulate Processes”*.

***Abstract:***
Wet granulation is a multiphase process utilised to produce aggregate particles with defined properties from very fine powders. Simulating this process on the microscale is challenging because of the large number of particles involved, which differ widely in both size and material properties. Macroscale methods, which track only the particle bulk properties, are efficient but do not resolve disperse particle properties such as the particle size distribution (PSD), which is key information for downstream processing. These deficiencies are addressed by mesoscale methods like population balance (PB) models, which track distributed properties such as the particle size by adding them as internal variables to the macroscale (CFD) model. However, most mesoscale methods are either inaccurate (method of moments when cutting off moments) or computationally expensive (Monte Carlo, class methods). Recently a new closure for the method of moments, the quadrature method of moments (QMOM), was introduced to allow accurate moment tracking of a PSD with low computational effort. Disadvantages of this method, e.g., it can suffer from instabilities,  can be overcome by non-dimensionalization. In this study we show our insights gained by non-dimensionalizing the QMOM equations for wet granulation processes, which model the PSD via growth, aggregation and breakage kernels. Relevant theoretical and numerical issues as well as limitations are discussed. Using constant kernels, the non-dimensionalised model is verified and validated thereafter on certain special cases. The effect of non-constant kernels on the non-dimensionalisation is discussed and we show that the non-dimensionalised model fails to accurately predict the moments of an experimental distribution using a volume-based solution to the PB equations. We show that a length-based solution to the PB equations can successfully predict the moments by setting non-constant kernels to fit the Sauter mean diameter.

***Description of the data in this software set:***
All data is provided in established open file formats (.csv, .md, .py and .svg). We refer to the paper for variable definition and names in generated 
figures from the software. Please do not move files or folders, as the data analysis tools are fully functional in the current state of this dataset.
The python files for the analysis are written in Python *3.7* and the MATLAB version used was *R2023a*. Python and MATLAB scripts are commented to make them self explanatory.


***Detailed description of the data and folder structure:***
Subfolders are named intuitively and the data inside subfolders should be described by the subfolders name.

```
Dataset/  
├── BuckinghamPy
│   ├── BuckinghamPy						: BuckinghamPy repository for version 1.0.3
│   └── BuckinghamTheorem.py					: Python script for the application of the Buckingham Theorem
└── QuadratureMethodOfMoments
    ├── Data 
    │   ├── InitialLactoseMCCPVP-CDF-v_L.csv			: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data
    │   └── InitialLactoseMCCPVP-n_V.csv			: Spreadsheet file with volume-based number particle size distribution (n_V) derived from InitialLactoseMCCPVP-CDF-v_L.csv
    │   ├── N13LactoseMCCPVP-CDF-v_L.csv 			: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N13
    │   └── N13LactoseMCCPVP-n_V.csv 	 			: Spreadsheet file with volume-based number particle size distribution (n_V) derived from N13LactoseMCCPVP-CDF-v_L.csv
    ├── Figures
    │   ├── AggregationAndBreakagedimensionlessM0M1.svg		: Figure showing the evolution of the dimensionless first and zeroth moment over time for aggregation and breakage
    │   ├── ApplicationMeanVolumeFit.svg			: Figure showing the fit of the mean particle volume for a twin-screw wet granulation experiment
    │   ├── ApplicationMomentComparison_nL.svg			: Figure showing logarithmic dimensionless length-based moments against moment order to compare the error qualitatively
    │   ├── ApplicationMomentComparison_nV.svg			: Figure showing logarithmic dimensionless moments against moment order to compare the error qualitatively
    │   ├── ApplicationSauterFit.svg				: Figure showing the fit of the Sauter mean diameter for a twin-screw wet granulation experiment
    │   ├── dimensionlessM0.svg					: Figure showing the evolution of the dimensionless zeroth moment over time
    │   ├── dimensionlessM1.svg					: Figure showing the evolution of the dimensionless first moment over time
    │   ├── PureAggregationdimensionlessM0M1.svg		: Figure showing the evolution of the dimensionless first and zeroth moment over time for pure aggregation
    │   ├── PureGrowthdimensionlessM0M1.svg			: Figure showing the evolution of the dimensionless first and zeroth moment over time for pure growth
    │   ├── M0.svg						: Figure showing the evolution of the zeroth moment over time
    │   └── M1.svg						: Figure showing the evolution of the first moment over time
    ├── Kernels
    │   ├── GrowthAggregationBreakage.m				: MATLAB function to compute the source terms for growth, aggregation and breakage
    │   └── GrowthHydrodynamicAggregationPowerLawBreakage.m	: MATLAB function to compute the source terms for growth, hydrodynamic aggregation and power law breakage
    ├── ApplicationQBMMPopulationBalanceNonDimensionalized.m	: MATLAB script for the application of the dimensionless Quadrature method of moment for a twin-screw wet granulation experiment using non-constant kernels
    ├── ComputeMoments.m					: MATLAB function to compute the moments from particle size distribution data
    ├── ConvertInitialDistribution.m				: MATLAB script which converts the initial distribution to a volume-based as well as length-based number distribution
    ├── QBMMPopulationBalance.m					: MATLAB script for the verification of the Quadrature method of moment for aggregation growth and breakage
    ├── QBMMPopulationBalanceNonDimensionalized.m		: MATLAB script for the verification of the dimensionless Quadrature method of moment for aggregation growth and breakage
    ├── ReadData.m						: MATLAB function that reads in a .csv file consisting of particle size distribution data
    ├──	 ValidationQBMMPopulationBalanceNonDimensionalized.m	: MATLAB script for the validation of the dimensionless Quadrature method of moment for special cases of aggregation, growth and breakage
    ├── Wheeler.m						: MATLAB function for the adaptive Wheeler algorithm
    ├── convertCDFtoPDF.m					: MATLAB function to convert a cumulative density function (CDF) into a probability density function (PDF)
    ├── convertPDFtoCDF.m					: MATLAB function to convert a probability density function (PDF) into a cumulative density function (CDF)
    ├── cutHighSizeRatio.m					: MATLAB function to cut high size ratios based on quantiles of the minimum and maximum particle sizes
    └── getMomenta.m						: MATLAB function to compute the moments from Gaussian quadrature weights and nodes
```