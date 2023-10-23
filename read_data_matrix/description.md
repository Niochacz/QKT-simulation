# Read data matrix
In this folder is posted code responsible for reading data matricec produced by other programs in project. It creates state vector map of probability plot or coherence from time function plot.

## Evolving function
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
Section where parameters are declared.

### Data matrix shape parameters
First you have to input what shape your data matrix is. Those parameters should be the same as those in file that created data matrix.

### Computation parameters
Parameters describing computation, which lead to given data matrix. Those are mainly used to read certain data file, while default name of saved file contains those parameters.
- N_Hlist - list of numbers of applied hamiltonians (applied sub-circuits)
- alphas - list of hamiltonian parameter
- thetalist2, philist2 - starting position of points
Those are all list, so program is capable of reading more than one file in one run.

### Coherence
Declaration of matrix to store coheretion from time relation.

## Execution section
Section, where user has to mark what program should produce from read data:
- plot - map of state vector probability
- coherence - plot of coherence from time relation

## Main function
Section resonsible for reading data, computing coherence and creating plots.  

First program iterate through list describing, which file to read. Then it computes classical evolution of starting point and load data to matrix. Last thing that program computes is so called coherence, which is integral over map of probability.

### Creating plot
If user marked to create plot then for each file in iteration a map of probability will be plotted and saved. For visualise reasons it is normalise by its integral.

### Coherence plot
If user marked to create coherence plot then for each *alpha* in  list *alphas* coherence function will be plotted and saved. Time in this function is described by number of applied hamiltonians.



