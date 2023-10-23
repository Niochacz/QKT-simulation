from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, IBMQ
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from qiskit.compiler import transpile


# Loading IBMQ account and setting backend
provider = IBMQ.load_account()
backend = provider.get_backend('ibmq_manila')

plot = False  # if creating phase space plot


# Evolving function
def evolve(theta, phi, p, k):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
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

# Declaration of list
circuits = []                                               # list with circuts object
jobs = []                                                   # list with jobs to do
resultlist = []                                             # list with probability outcomes


# Circuit creation

# Setting number of qubits and bits per circuit
qreg = QuantumRegister(3, 'q')
creg = ClassicalRegister(3, 'c')


# Sub-circuit as Hamiltonian implementation
alpha = 4*np.pi/8                                           # hamiltonian parameter
kick = QuantumCircuit(qreg, name = 'kick')                  # declaration of circuit

# Appending gates
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


# Setting point starting position
theta0, phi0 = 2*np.pi/5, 3*np.pi/4


for phi in philist:
    for theta in thetalist:
        qc = QuantumCircuit(qreg, creg)                     # main circuit creation
        
        # rotating qubits to starting position
        qc.ry(-theta0, qreg)
        qc.rz(-phi0, qreg)
        
        # Appending sub-circuits to main circuit
        for i in range(N_H):
            qc.append(kick, qreg)
            
        # rotating qubits to |theta, phi> state
        qc.rz(phi, qreg)
        qc.ry(theta, qreg)
        
        # probability measure
        qc.measure(qreg, creg)
        
        # Appending main circuit to list
        circuits.append(qc)
        
        
Ncirc = len(thetalist)*len(philist)                         # number of all circuits (one per pixel)


# Circuit execution

# Running circuits on IBM computers (take a while to wait in queue [I mean ~1 hour is not that bad])
for i in range(int(np.floor(Ncirc/100))):
    jobs.append(backend.run(transpile(circuits[i*100:(i+1)*100], backend), shots = shots))
LastJob = backend.run(transpile(circuits[int(np.floor(Ncirc/100))*100:], backend), shots = shots)


# Retrieving number of counts to matrix
for i, job in enumerate(jobs):
    for circuit in circuits[i*100:(i+1)*100]:
        counts = job.result().get_counts(circuit)
        if '000' in counts.keys():
            resultlist.append(counts['000']/shots)
        else:
            resultlist.append(0)

for circuit in circuits[int(np.floor(Ncirc/100))*100:]:
    counts = LastJob.result().get_counts(circuit)
    if '000' in counts.keys():
        resultlist.append(counts['000']/shots)
    else:
        resultlist.append(0)


# Reshaping result list to plotable format
resultlist = np.asarray(resultlist).reshape(len(philist), len(thetalist))
resultlist = np.transpose(resultlist)

# Saving data matrix
np.savetxt(f'man_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.txt', resultlist)


# Creating plot
if plot:
    
    # Computing classical evolution
    theta1, phi1 = theta0, phi0                                # iterator
    for i in range(N_H):
        theta1, phi1 = evolve(theta1, phi1, np.pi/2, -alpha)
        
    plt.figure(figsize=(4, 2))
    ax = plt.gca()
    im = ax.imshow(resultlist, vmin = 0.0, vmax = 1.0, cmap = 'inferno')
    plt.plot((16+phi0/phistep) % 31, (14-theta0/thetastep) % 15, 'wo')
    plt.plot((16+phi1/phistep) % 31, (14-theta1/thetastep) % 15, 'go')
    plt.xlabel('$\phi$')
    plt.xticks([0, int(len(philist)/2), len(philist)-1], ['$-\pi$', '0', '$\pi$'])
    plt.ylabel('$\Theta$')
    plt.yticks([0, len(thetalist)-1], ['$\pi$', '0'])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax = cax)
    plt.tight_layout()
    plt.savefig(f'man_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.png')
