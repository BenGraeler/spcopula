# spcopula

The spcopula package allows to model spatial and spatiotemporal random fields. 
The dependence structure is captured with the help of Vine Copulas in which the bivariate building blocks depend on the separating distances between the locations. 
The change in the strength of dependence can be continuously modelled through a correlogram (also compare variogram in geostatistics/kriging) based on Kendall's tau. 
The functionality includes the estimation of the dependence structure and the interpolation and simulation of the modelled random fields. 
Additionally, the spcopula package features functions to calculate multivariate return periods based on bivariate copulas or vine copulas.

The package has initially been developed through the DFG project "Developing Spatio-Temporal Copulas" at the Institute for Geoinformatics, University of Münster.
Its state can be rated as a "proof of concept", as it has evolved along a series of use-cases. 
However, it is supposed to be applicable to many other use-cases beyond these and has been used in a number of scientific publications.
