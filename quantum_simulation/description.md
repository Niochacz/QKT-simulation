# Quantum Kicked Top simulation on IBM classical computer
In this folder is posted code responsible for collecting data about act of Quantum Kicked Top hamiltonian on 3 qubits computed on IBM classical computer.

## Setting backend
In this program we will use 'qasm_simulator' backend. It is backend provided by IBM that simulate quantum computer on classical computer.

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
First you have to set what size your outcome data matrix should be. The lesser step, the bigger plot resolution, but also greater computation time.

### Computation parameters
There are 5 parameters:
- shots - number of sample per pixel (how many times one circuit will be executed to best estimate probability)
- N_H - number of applied hamiltonian on qubits
- thetalist2, philist2 - list of starting points position
- alphas - list of hamiltonian parameter

## Executing section
Section, where user has to mark what program should produce from computed data:

- plot - map of state vector probability
- savetotxt - save outcome data matrix to txt file

## Main function
Program can produce in one run a couple of data. The price for such comfort is giant loop, because everything from this point needs to be simulated for certain alpha. Fortunately simulations does not take that long, so it is sufficient solution.

## Circuit creation

### Setting qubit and bit number
Almost last parameter to set in our experiments is number of qubits.
It depends on what are you trying to measure or simulate.
In our experiment that number is set to 3 and hamiltonian implementation is also built for 3 qubits.

### Creating hamiltonian implementation
Next there is implemented QKT Hamiltonian for 3 qubits as sub-circuit.
First, hamiltonian parameter and such object as circuit is declared.
Then to declared circuit gates are being attached. Such circuit has explicit action as QKT Hamiltonian.

### List declaration
- circuits - list of all to be executed circuits
- jobs - list of all jobs (circuit execution)
- resultlist - list with number of successful shots for each job

### Main circuit creation  
Last parameter to declare in our experiment is starting position of evolving point (position in spherical coordinate system).  
Then for each pixel on plot there is creating a circuit, which consist of:  
- rotating qubits to starting position
- appending sub-circuits in number as declared in N_H (number of Hamiltonians)
- rotating qubits to |theta, phi> state
- measurement

## Circuit execution
In this section circuits are uploaded to compute by IBM computers. Once they are executed number of counts (successful shots) for each circuit is being retrieved into a matrix and normalised to be treated as probability. 

## Save to txt file
If mark as to be done then last section saves the outcome data matrix to .txt file.

## Plot creation
If mark as to be done then last section produces plot from computed data.
Plotted is plane corresponding to "Bloch sphere" for 3 qubits. It is a map of probability, where state vector will be.
Also on graph is plotted starting satte vector position (white dot) and classical evolution of starting point (green dot).




