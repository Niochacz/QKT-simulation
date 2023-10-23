from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, Aer, IBMQ
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from qiskit.compiler import transpile
from qiskit.providers.aer.noise import NoiseModel


# Loading IBMQ account and setting backend
provider = IBMQ.load_account()
backend = Aer.get_backend('qasm_simulator')
noise_model = NoiseModel.from_backend(provider.get_backend('ibmq_manila'))
coupling_map = backend.configuration().coupling_map
basis_gates = noise_model.basis_gates


# Evolving function
def evolve(theta, phi, p, k):
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)

    dx = (x*np.cos(p) + z*np.sin(p)) * np.cos(k*(z*np.cos(p) - x*np.sin(p))) - y*np.sin(k*(z*np.cos(p) - x*np.sin(p)))
    dy = (x*np.cos(p) + z*np.sin(p)) * np.sin(k*(z*np.cos(p) - x*np.sin(p))) + y*np.cos(k*(z*np.cos(p) - x*np.sin(p)))
    dz = -x*np.sin(p) + z*np.cos(p)

    dtheta = np.arccos(dz)
    dphi = np.arctan2(dy, dx)
    return dtheta, dphi


# Initial parameters

# First 4 parameters are describing outcome shape (data matrix shape).
# Better do not to change or set same in other files in project
thetastep = np.pi/15
phistep = np.pi/16
thetalist = np.flipud(np.arange(0, np.pi, thetastep))       # list with angles theta from pi to 0
philist = np.arange(-np.pi, np.pi, phistep)                 # list with angles phi from -pi to pi

# Computation parameters
shots = 10000                                               # number of shots per one pixel on plot
N_H = 9                                                     # number of kicks (applied hamiltonians)
thetalist2 = [2/5*np.pi]                                    # list of starting point theta position
philist2 = [3/4*np.pi]                                      # list of starting point phi position
alphas = [1/8*np.pi, 2/8*np.pi, 4/8*np.pi]                  # list of hamiltonian parameter


# Executing sections
savetotxt = True                                            # if save outcome data matrix to txt file
plot = False                                                # if create a map of probability plot


# Main function
for alpha in alphas:
    
    # Circuit creation
    
    # Setting number of qubits and bits per circuit
    qreg = QuantumRegister(3, 'q')
    creg = ClassicalRegister(3, 'c')

    # Sub-circuit as Hamiltonian implementation
    kick = QuantumCircuit(qreg, name = 'kick')
    kick.ry(-np.pi/2, qreg)
    kick.cnot(qreg[1], qreg[0])
    kick.rz(alpha, qreg[0])
    kick.cnot(qreg[2], qreg[1])
    kick.rz(alpha, qreg[1])
    kick.cnot(qreg[0], qreg[1])
    kick.rz(alpha, qreg[1])
    kick.cnot(qreg[2], qreg[1])
    kick.cnot(qreg[0], qreg[1])
    kick.cnot(qreg[1], qreg[0])

    for phi0 in philist2:
        for theta0 in thetalist2:
            
            # Declaration of list
            circuits = []                                               # list with circuts object
            jobs = []                                                   # list with jobs to do
            resultlist = []                                             # list with probability outcomes

            for phi in philist:
                for theta in thetalist:
                    # Main circuit creation
                    qc = QuantumCircuit(qreg, creg) 
                    
                    # rotating qubits to starting position
                    qc.ry(-theta0, qreg)
                    qc.rz(-phi0, qreg)

                    # appending sub-circuits to main circuit
                    for i in range(N_H):
                        qc.append(kick, qreg)
                    
                    # rotating qubits to |theta, phi> state
                    qc.rz(phi, qreg)
                    qc.ry(theta, qreg)
                    
                    # measurement
                    qc.barrier(qreg)
                    qc.measure(qreg, creg)
                    
                    circuits.append(qc)
            Ncirc = len(thetalist)*len(philist)

            # Circuit execution
            
            # Running circuits on IBM classical computers
            jobs.append(backend.run(transpile(circuits, backend), shots = shots, coupling_map=coupling_map, basis_gates=basis_gates, noise_model=noise_model))

            # Retrieving number of counts to matrix
            for i, job in enumerate(jobs):
                for circuit in circuits:
                    counts = job.result().get_counts(circuit)
                    if '000' in counts.keys():
                        resultlist.append(counts['000']/shots)
                    else:
                        resultlist.append(0)
                        
            # Reshaping result list to plotable format
            resultlist = np.asarray(resultlist).reshape(len(philist), len(thetalist))
            resultlist = np.transpose(resultlist)
            
                
            # Saving to txt file
            if savetotxt:
                np.savetxt(f'sym_noise_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.txt', resultlist)
            
                
            # Creating plot
            if plot:
                
                # Computing classical evolution
                thetan, phin = theta0, phi0
                for i in range(N_H):
                    thetan, phin = evolve(thetan, phin, np.pi/2, -alpha)
                    
                plt.figure(figsize=(4, 2))
                ax = plt.gca()
                im = ax.imshow(resultlist, vmin = 0.0, vmax = 1.0, cmap = 'inferno')
                plt.xlabel('$\phi$')
                plt.xticks([0, int(len(philist)/2), len(philist)-1], ['$-\pi$', '0', '$\pi$'])
                plt.ylabel('$\Theta$')
                plt.yticks([0, len(thetalist)-1], ['$\pi$', '0'])
                ax.plot((16+phi0/phistep) % 31, (14-theta0/thetastep) % 15, 'wo')
                ax.plot((16+phin/phistep) % 31, (14-thetan/thetastep) % 15, 'go')
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1)
                plt.colorbar(im, cax = cax)
                plt.tight_layout()
                plt.savefig(f'sym_noise_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.png')
                plt.clf()
