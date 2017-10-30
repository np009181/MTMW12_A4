#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
import numpy as np

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
    
    #setting new parameters (incl. derived) to use for q5#
    d = 0.32
    dx = (xmax - xmin)/(nx - 1)    
    dt = d*dx**2/K
         #Defining the non-dim diffusion coefficient, K*dt/dx**2 but
               #now we fix value for q5#
               
    #the following code allows us to perform code with same time duration#

    NT = int(4/dt)
    if (4/dt - NT < 0.5):
        nt = NT + 1
    else:
        nt = NT
    
    #New values for K so values of d are changed, for q4#
    K1 = 3e-3
    K2 = 3.5e-3
    K3 = 2.5e-3
         
    
    d1 = K1*dt/dx**2
    d2 = K2*dt/dx**2
    d3 = K3*dt/dx**2
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
    
    #Finding the diffusion using FTCS and BTCS with different nt#
    phiFTCS1 = FTCS(phiOld.copy(), d, nt = 400)
    phiBTCS1 = BTCS(phiOld.copy(), d, nt = 400)
    
    #Finding diffusion using FTCS and BTCS with new value of d#
    phiFTCS2 = FTCS(phiOld.copy(), d1, nt)
    phiBTCS2 = BTCS(phiOld.copy(), d1, nt)
    phiFTCS3 = FTCS(phiOld.copy(), d2, nt)
    phiBTCS3 = BTCS(phiOld.copy(), d2, nt)
    phiFTCS4 = FTCS(phiOld.copy(), d3, nt)
    phiBTCS4 = BTCS(phiOld.copy(), d3, nt)
    
    #Working out the error norms and printing them#
    
    print("FTCS L2 error norm = ", L2ErrorNorm (phiFTCS, phiAnalytic))
    print("BTCS L2 error norm =", L2ErrorNorm (phiBTCS, phiAnalytic))
    
    #Now we plot the solutions#
    font = {'size'  : 15}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='blue')
    plt.plot(x, phiBTCS, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/FTCS_BTCS.pdf')
    
    #Finding the error of the numerical solutions#
    eFTCS = phiAnalytic - phiFTCS
    eBTCS = phiAnalytic - phiBTCS
    
    font = {'size'  : 10}
    plt.rc('font', **font)
    plt.figure(2)
    plt.plot(x, eFTCS, label='Error FTCS', color='blue')
    plt.plot(x, eBTCS, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/Error_FTCS_BTCS.pdf')
    
    font = {'size'  : 15}
    plt.rc('font', **font)
    plt.figure(3)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS1, label='FTCS', color='blue')
    plt.plot(x, phiBTCS1, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/New_FTCS_BTCS.pdf')
    
    #Working out the error for FTCS1 and BTCS1#
    eFTCS1 = phiAnalytic - phiFTCS1
    eBTCS1 = phiAnalytic - phiBTCS1
    
    font = {'size'  : 10}
    plt.rc('font', **font)
    plt.figure(4)
    plt.plot(x, eFTCS1, label='Error FTCS', color='blue')
    plt.plot(x, eBTCS1, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/NewError_FTCS_BTCS.pdf')
    
    font = {'size'  : 15}
    plt.rc('font', **font)
    plt.figure(5)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS2, label='FTCS', color='blue')
    plt.plot(x, phiBTCS2, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/ChangeD_FTCS_BTCS.pdf')
    
    font = {'size'  : 15}
    plt.rc('font', **font)
    plt.figure(6)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS3, label='FTCS', color='blue')
    plt.plot(x, phiBTCS3, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/ChangeD1_FTCS_BTCS.pdf')
    
    #Working out the error of FTCS and BTCS with new d value#
    eFTCS2 = phiAnalytic - phiFTCS2
    eBTCS2 = phiAnalytic - phiBTCS2
    
    eFTCS3 = phiAnalytic - phiFTCS3
    eBTCS3 = phiAnalytic - phiBTCS3
    
    eFTCS4 = phiAnalytic - phiFTCS4
    eBTCS4 = phiAnalytic - phiBTCS4
    
    font = {'size'  : 10}
    plt.rc('font', **font)
    plt.figure(7)
    plt.plot(x, eFTCS2, label='Error FTCS', color='blue')
    plt.plot(x, eBTCS2, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/D_Error_FTCS_BTCS.pdf')
    
    font = {'size'  : 10}
    plt.rc('font', **font)
    plt.figure(8)
    plt.plot(x, eFTCS3, label='Error FTCS', color='blue')
    plt.plot(x, eBTCS3, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.1,1))
    plt.xlabel('$x$')
    plt.savefig('plots/D1_Error_FTCS_BTCS.pdf')
    
def OrderofConvergence (N, fixD):
    
    #Code to initialise the vector to store the values for L2 error norm#
    l2 = np.zeros((N, 3))
    
    for it in range (N):
        l2[it,:] = main(nx=21+it*50, d=fixD)
        
    # Introducing F and B as FTCS and BTCS gradients between individual#
    # points respectively. This should give an vector for both#
    F = np.zeros(N - 1)
    for i in range(N - 1):
        F[i] = (np.log(l2[i + 1, 1]) - np.log(l2[i, 1]))/(np.log(l2[i + 1, 0]) - np.log(l2[i,0]))
    
        
    B = np.zeros(N - 1)
    for i in range(N - 1):
        B[i] = (np.log(l2[i + 1, 2]) - np.log(l2[i, 2]))/(np.log(l2[i + 1, 0]) - np.log(l2[i,0]))
    
    return OrderofConvergence

    #Graph to show the Log-Log graph#
    
    font = {'size'  : 10}
    plt.rc('font', **font)
    plt.clf()
    plt.ion()
    plt.figure(9)
    plt.loglog(l2[:,0].transpose(), l2[:,1].transpose,
               label='FTCS (mean slope = {:.2f})'.format(n.mean()))
    plt.loglog(l2[:,0].transpose(), l2[:,2].transpose(),
               label='BTCS (mean.slope = {:.2f})'.format(m.mean()))
    plt.loglog(l2[:,0].transpose(),8*(l2[:,0].transpose())**2)
    plt.legend()
    plt.xlabel('Delta x')
    plt.ylabel('L2 norm error')
    plt.title('Log-Log plot showing L2 norm error against a fixed D')
    plt.savefig('plots/L2_Error{}.pdf')
    
    
    
main()
