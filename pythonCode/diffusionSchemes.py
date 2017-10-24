# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np

# The linear algebra package for BTCS (for solving the matrix equation)
import scipy.linalg as la

def FTCS(phiOld, d, nt):
    nx = len(phiOld)
    
    #Create a new time-step array for phi#
    phi = phiOld.copy()
    

   #FTCS for all time steps, need to add formula for endpoints#
   for it in xrange(int(nt)):
        phi[0] = phiOld[0] + d*(phiOld[2] - 2*phiOld[1] + phiOld[0])
        
        for j in range (1, nx - 1):
            phi[j] = phiOld[j] + d*(phiOld[j] - 2*phiOld[j] + phiOld[j])
            phi[nx - 1] = phiOld[0] + d*(phiOld[nx - 3] - 2*phiOld[nx - 2] + phiOld[nx - 1])
            
        phiOld = phi.copy() 
            

    return phi
        


def BTCS(phi , d , nt) :
    nx = len (phi)
    
# array representing BTCS 
    M = np . zeros ([ nx , nx ]) 
    
# Zero gradient boundary conditions 
    M[0 ,0] = 1.
    M[0 ,1] = -1.
    M[-1,-1] = 1.
    M[-1,-2] = -1.
    
    
    for i in xrange (1 , nx - 1):
        
        M[i , i -1] = -d
        M[i , i ] = 1+2*d
        M[i , i +1] = -d
# BTCS for a l l time steps
    for it in xrange (int(nt)):
     
     # RHS for zero gradient boundary conditions 
        phi [0] = 0
        phi [-1] = 0
     
     
        phi = la . solve (M, phi )


    return phi
