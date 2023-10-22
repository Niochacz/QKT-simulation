# Classical Kicked Top evolution  
In this folder is posted code responsible for plotting phase space of classical Kicked Top hamiltonian.

## Evolving funtion  
This is implementation of classical evolution mainly based on this paper (Section III):   
> Udaysinh T. Bhosale and M. S. Santhanam, Periodicity of quantum correlations in the quantum kicked top, Phys. Rev. E 98, 052228 – Published 28 November 2018. DOI:https://doi.org/10.1103/PhysRevE.98.052228

But also simlarly to archetype of QKT (Y and Z axis are swapped) [Section 3]:  
> Kuś, M.; Scharf, R. and Haake, F (1987). Symmetry versus degree of level repulsion for kicked quantum systems. Z. Phys. B 66: 129. doi:10.1007/bf01312770   

Function takes arguments about current position of point and returns position of given point in next step. 
 - takes position of point in spherical coordinate system
 - transform to cartesian coordinates
 - computes position of point in next step following equations in mentioned papers
 - returns position in spherical coordinate system

*p* and *k* are hamiltonian parameters.

## Initial parameters  
Section where initial parameters are defined:
- steps - number of evolution steps
- p0, k0 - hamiltonian parameters
- Ntheta - number of grid points along theta axis on plot
- Nphi - number of grid points along phi axis on plot
- thetas, phis - matrices of computed positions in each step

## Initial start points  
Creates a grid of evenly spaced points.  

## Main function  
Section where position of each point is being computed.  

## Creating plot  
Section responsible for creating plot of phase space.
The outcome is a scatter plot of all points positions across evolution. 
