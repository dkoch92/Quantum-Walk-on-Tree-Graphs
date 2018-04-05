import Ntreefunctions as nt
import matplotlib.pyplot as plt
import U_Eigenvalues as egn
import Regression_Fits as rg
from numpy import *

'''
This code is for generating A(M),B(M),Rho,Gamma for the function epsilon(N,M)

The user specifies the min/max values for N and M

The code generates the matrix U's, finds their eigenvalues/eigenvectors via U_Eigenvalues functions, 
caluclates the leading order eigenangles, and uses them in best fits

The code calls on the module Regression_Fits.py for fitting / plotting
Final constants for epsilon(N,M) are printed
'''

#-------------- Setting the Size --------------------------------
N_max = 14     
N_min = 2     #Define the min / max graph sizes (not recommended going over N=10 / M=15)
M_max = 8     #Minimum N/M that will work is 2
M_min = 2

#Creating Arrays / Matrices to store values
Angles   = zeros(shape=(N_max-N_min+1,M_max-M_min+1))
Betas    = zeros(shape=(N_max-N_min+1,M_max-M_min+1))
N_axis   = zeros(N_max-N_min+1)                          #arrays for plotting
M_axis   = zeros(M_max-M_min+1)                

#-------------- Running The Code --------------------------------


for i in arange(N_min,N_max+1):
    N = int(i)
    N_axis[i-N_min] = N
    for j in arange(M_min,M_max+1):
        #-----------Generating U Matrices--------------
        M = int(j)  
        M_axis[j-M_min] = M    
        Int,eigen = egn.Find_Eigen(N,M)
        #-------------Extracting Info-----------------
        B,B_r,B_angle = egn.Overlaps(Int,eigen)                             #stores the overlaps of each eiegenstate with the initial state of the system
        Estates_r,Estates_angle = egn.Eigenstates_Polar(eigen[1])           #stores all the eiegenstates as polar values
        Evalues_angle           = egn.Eigenvalues_Polar(eigen[0])           #stores all the eigenvalues as polar values
        top2 = egn.Find_Largest(B_r)                                        #finds the lcoation of the largest two Beta's
        index = int(top2[0])
        #------------                              
        Angles[i-N_min,j-M_min]   = abs(Evalues_angle[index])
        Betas[i-N_min,j-M_min]    = B_r[index]

#----------- Power Fits for Constant M ----------------------------
Am = zeros( len(M_axis))
Bm = zeros( len(M_axis))
rm = zeros( len(M_axis))

plot_choices = array([2,3,4,8])             #select which values of M to plot

for i in arange( len( M_axis) ):
    Am[i],Bm[i],rm[i] = rg.Power_Regression(N_axis,Angles[:,i])     #calculates the A,B constants as well as r (correlation coefficient)



#-----------------Best Fits for A(M) and B(M)---------------------
alpha,beta,ra = rg.Power_Regression(M_axis,Am)
gamma,rho,rb  = rg.Linear_Regression(M_axis,Bm)



#-----------------Plotting The Results ---------------------------
rg.Power_Fit_Plot(M_axis,N_axis,Angles,Am,Bm,rm,plot_choices)    #comment this out out if you do not wish to plot   
print("_____Constants for epsilon(N,M) (in degrees)_____")
print("Alpha:  ",round(alpha,5),"     Beta:  ",round(beta,5))
print("Gamma:  ",round(gamma,5),"      Rho:  ",round(rho,5))


