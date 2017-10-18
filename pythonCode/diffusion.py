#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
execfile("diffusionSchemes.py")
execfile("diagnostics.py")
execfile("initialConditions.py")

def main():
    #Parameters
    xmin = 0
    xmax = 1
    nx = 41
    nt = 40
    dt = 0.1
    K = 2e-3
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    
    #Derived parameters#
    dx = (xmax - xmin)/(nx - 1) #We want 40 time-steps so we have to -1#    
    d = K*dt/dx**2   #Defining the non-dim diffusion coefficient#
    
    #Working out the spacial points#
    x = np.zeros(nx)
    for j in xrange (nx):
        x[j] = xmin + j*dx
    print ('x=', x)

    #Initial conditions#
    phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    #Analytic solution of square profile in inf domain#
    phiAnalytic = analyticErf(x, K*dt*nt, squareWaveMin, squareWaveMax)
    
    #Finding diffusion using FTCS and BTCS#
    phiFTCS = FTCS(phiOld.copy(), d, nt)
    phiBTCS = BTCS(phiOld.copy(), d, nt)
    
    #Working out the error norms and printing them#
    print("FTCS L2 error norm = ", L2ErrorNorm (phiFTCS, phiAnalytic))
    print("BTCS L2 error norm =", L2ErrorNorm (phiBTCS, phiAnalytic))
    
    #Now we plot the solutions#
    font = {'size'  : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='__', linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='blue')
    plt.plot(x, phiBTCS, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/FTCS_BTCS.pdf')
    
    
    
main()
    
    
