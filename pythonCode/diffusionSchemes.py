# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np

# The linear algebra package for BTCS (for solving the matrix equation)
import scipy.linalg as la

def FTCS(phiOld, d, nt):
    nx = len(phiOld)
    
    #Create a new time-step array for phi#
    phi = phiOld.copy()
    
    #FTCS for all time steps#
    for it in xrange(int(nt)):
        phi[0] = ???
        
        
        
    for i in range (1, nx - 1):
        phi[i] = ???
