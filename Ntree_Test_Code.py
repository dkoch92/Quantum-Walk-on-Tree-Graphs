#------------- Modules from Github Project ---------------
import Ntreefunctions as nt
import Quantum_Walk_Simulation as qw
import Regression_Fits as rg
import U_Eigenvalues as egn
#-------------Standard Python Libraries ------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import linalg as LA
import scipy as sci
import itertools as it
'''
This code is meant to test that all python modules load correctly

Nothing in this code shuld be altered, simply run the code

Make sure that this .py code is running in the same directory as all the other modules:
Ntreefunctions.py | Quantum_Walk_Simulation.py | Regression_Fits.py | U_Eigenvalues.py


'''
#-----------Checking Everything----------------

if( nt and qw and rg and egn ):
    print(" All N-Tree Modules found ")         #if this messege prints, all the N-Tree modules from github are correctly located
else:
    print(" Not all N-Tree Modules found ")     #you will most likely get an error before reachign this line
        
    
if( plt and Axes3D and np and LA and sci and it ):
    print(" All standard Python libraries found ")      #if this messege prints, all standard python libararies loaded correctly
else:
    print(" Not Python libraries found ")               #you will most likely get an error before reaching this line
        