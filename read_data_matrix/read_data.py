import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
Program to read data produced by quantum_computation.py, quantum_simulation.py
"""


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

# First 4 parameters are describing data to read matrix, therefore they should be equal to those used to create data file.
# Better do not to change or set same in other files in project
thetastep = np.pi/15
phistep = np.pi/16
thetalist = np.flipud(np.arange(0, np.pi, thetastep))     # list with angles theta from pi to 0
philist = np.arange(-np.pi, np.pi, phistep)               # list with angles phi from -pi to pi

# Those parameters are used to determine which file should be read
N_Hlist = range(11)                                       # number of kicks (applied hamiltonians)
alphas = [1/8*np.pi, 2/8*np.pi, 4/8*np.pi]                # hamiltonian parameter
thetalist2 = [2/5*np.pi]                                  # theta starting position of points
philist2 = [3/4*np.pi]                                    # phi starting position of points

# Last parameter is used to compute coherence
totQlist = np.zeros((len(alphas), len(N_Hlist)))          # coherence value matrix


# Executing sections
plot = True                                               # if create a phase space plot
coherence = False                                         # if create a coherence plot


# Main function
for a, alpha in enumerate(alphas):
    for N, N_H in enumerate(N_Hlist):
        totQ = 0
        for phi0 in philist2:
            for theta0 in thetalist2:
                thetan, phin = theta0, phi0
                for i in range(N_H):
                    thetan, phin = evolve(thetan, phin, np.pi/2, -alpha)  # classicly evolved point  
                
                # Loading data
                resultlist = np.loadtxt(f'./man_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.txt')
                
                # computing coherence
                for i, theta in enumerate(thetalist):
                    totQ += np.sum(resultlist[i])*np.sin(theta)
                totQlist[a, N] = totQ*(2*np.pi/(len(philist)*len(thetalist)))
                
                if plot:  # create plot or not
                
                    resultlist *= 1/totQlist[a, N]  # normalising
                    
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
                    cax = divider.append_axes("right", size="5%", pad = 0.1)
                    plt.colorbar(im, cax = cax)
                    plt.tight_layout()
                    plt.savefig(f'man_domnozone_{np.round(alpha,2)}_{np.round(theta0,2)}x{np.round(phi0,2)}_{N_H}.png')
                    plt.clf()

if coherence:  # to create coherence plot or not
    plt.plot(N_Hlist, totQlist[0], 'bo')
    plt.plot(N_Hlist, totQlist[1], 'go')
    plt.plot(N_Hlist, totQlist[2], 'ro')
    plt.xlabel('Number of kicks $M$')
    plt.ylabel('Value of $I$')
    plt.legend([r'$\alpha = \frac{\pi}{8}$', r'$\alpha = \frac{\pi}{4}$', r'$\alpha = \frac{\pi}{2}$'])
    plt.show()
