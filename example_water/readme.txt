###
Eg_water.R: Case study: H2O potential energy (without PCA).

Eg_water_pca.R: Case study: H2O potential energy (with PCA).

chemrev_nuprime-theta-grid_computed.xyz: Contain H2O molecular configurations arranged on a grid, with varying angles and asymmetric stretches. Dataset from the manuscript "Gaussian Process Regression for Materials and Molecules" (submitted to Chemical Reviews).

pswater-ipi.out: Contain energies (including the dipole) computed with the Partridge–Schwenke model by i-PI. Dataset from the manuscript "Gaussian Process Regression for Materials and Molecules" (submitted to Chemical Reviews).

water.ipynb: Jupyter notebook with the code for Section 3.3, "H2O Potential Energy: A Hands-On Example," from the manuscript "Gaussian Process Regression for Materials and Molecules" (submitted to Chemical Reviews).

water2.ipynb: Jupyter notebook that uses the datasets chemrev_nuprime-theta-grid_computed.xyz and pswater-ipi.out to construct the feature matrix X_soap.csv and the corresponding energy file y_energy.csv through Smooth Overlap of Atomic Positions (SOAP) method.

X_soap.csv: SOAP feature matrix for all configurations.

y_energy.csv: SOAP feature matrix for all configurations.
###
