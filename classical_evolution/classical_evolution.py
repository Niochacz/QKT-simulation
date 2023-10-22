import numpy as np
import matplotlib.pyplot as plt


#Evolving function
def evolve(theta, phi, p, k):
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)

    dx = (x*np.cos(p) + z*np.sin(p)) * np.cos(k*(z*np.cos(p) - x*np.sin(p))) - y*np.sin(k*(z*np.cos(p) - x*np.sin(p)))
    dy = (x*np.cos(p) + z*np.sin(p)) * np.sin(k*(z*np.cos(p) - x*np.sin(p))) + y*np.cos(k*(z*np.cos(p) - x*np.sin(p)))
    dz = -x*np.sin(p) + z*np.cos(p)

    dtheta = np.arccos(dz)
    dphi = np.arctan2(dy,dx)
    return dtheta, dphi


#Initial parameters
steps = 100
p0 = np.pi/2
k0 = 3*np.pi/8
Ntheta = 30
Nphi = 60
thetas = np.zeros((Ntheta, Nphi, steps))
phis = np.zeros((Ntheta, Nphi, steps))


#Initial start points
for i in range(Ntheta):
    thetas[i,:,0] = [i*np.pi/Ntheta for p in range(Nphi)]
    for j in range(Nphi):
        phis[:,j,0] = [-np.pi + 2*j*np.pi/Nphi for p in range(Ntheta)]
        

#Main function      
for i in range(Ntheta):   
    for j in range(Nphi): 
        theta0 = thetas[i,j,0]
        phi0 = phis[i,j,0]
        for n in range(1, steps):
            theta0, phi0 = evolve(theta0, phi0, p0, k0)
            thetas[i,j,n] = theta0
            phis[i,j,n] = phi0


#Creating plot
phis = np.reshape(phis, (Ntheta*Nphi*steps))
thetas = np.reshape(thetas, (Ntheta*Nphi*steps))
plt.plot(phis,thetas, 'b.', ms = 0.1)
plt.xlabel('$\phi$')
plt.xlim(-np.pi, np.pi)
plt.ylabel('$\Theta$')
plt.ylim(0,np.pi)
plt.title(f'Phase space for p = {round(p0,2)} and k = {round(k0,2)}')
plt.show()
