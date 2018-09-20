#! /usr/bin/python3
import numpy as np
from scipy import linalg as la 
import matplotlib
import matplotlib.pyplot as plt

############################
print("Loading matrices...")
############################
A  = np.loadtxt('A.txt')
nc = int(np.sqrt(A.shape[0]))
A  = np.reshape(A,(nc,nc))
P  = np.loadtxt('P.txt')
nc = int(np.sqrt(P.shape[0]))
P  = np.reshape(P,(nc,nc))

###################################
print("Computing eigenvalues...\n")
###################################
val = la.eigh(A,P,eigvals_only=True)
absval = np.sort(np.abs(val))
np.savetxt('eigvals.txt',absval,delimiter='\n')
plt.plot(absval,'-o')
plt.grid()
plt.show()


















