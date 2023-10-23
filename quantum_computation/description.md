# Quantum Kicked Top simulation on IBM quantum computer
In this folder is posted code responsible for collecting data about act of Quantum Kicked Top hamiltonian on 3 qubits computed on IBM quantum computer.

## Loading IBMQ account
To use IMBQ computers you have to create your account. Every account has assigned token, which you should provide to use this code.
>https://qiskit.org/ecosystem/ibmq-provider/tutorials/1_the_ibmq_account.html
Next you have to select backend, which you will be using.

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
The resolution of plot is also limited by numbers of circuit you can upload to queue with free IBMQ account. The limit is 500 circuits, so by default one of the highest resolution is set (32,15).

### Computation parameters
There are 2 parameters:
- shots - number of sample per pixel (how many times one circuit will be executed to best estimate probability)
- N_H - number of applied hamiltonian on qubits

### List declaration
- circuits - list of all to be executed circuits
- jobs - list of all jobs (circuit execution)
- resultlist - list with number of successful shots for each job

## Circuit creation

### Setting qubit and bit number
Almost last parameter to set in our experiments is number of qubits. 
It depends on what are you trying to measure or simulate.
In our experiment that number is set to 3 and hamiltonian implementation is also built for 3 qubits.
Be aware also that you have to obey IBMQ computers architecture.

### Creating hamiltonian implementation
Next there is implemented QKT Hamiltonian for 3 qubits as sub-circuit.  
First, hamiltonian parameter and such object as circuit is declared.  
Then to declared circuit gates are being attached. Such circuit has explicit action as QKT Hamiltonian.

### Main circuit creation  
Last parameter to declare in our experiment is starting position of evolving point (position in spherical coordinate system).  
Then for each pixel on plot there is creating a circuit, which consist of:  
- rotating qubits to starting position
- appending sub-circuits in number as declared in N_H (number of Hamiltonians)
- rotating qubits to |theta, phi> state
- measurement

## Circuit execution
In this section circuits are uploaded to queue in IBMQ servers. Once they are executed number of counts (successful shots) for each circuit is being retrieved into a matrix and normalised to be treated as probability. In the end the outcome data matrix is saved to .txt file.
(LastJob is out of loop, becouse of some sort of bug, so the last iteration is done manually)

## Plot creation
If mark as to be done then last section produces plot from computed data.
Plotted is plane corresponding to "Bloch sphere" for 3 qubits. It is a map of probability, where state vector will be.
Also on graph is plotted starting satte vector position (white dot) and classical evolution of starting point (green dot).

