import Ntreefunctions as nt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import U_Eigenvalues as egn
from numpy import *

'''
This code is used for veiwing the function epsilon(N,M) as a 3D plot

The user specifies the min/max values for N and M

The code generates the matrix U's, finds their eigenvalues/eigenvectors via U_Eigenvalues functions, 
caluclates the leading order eigenangles, and plots them as a 3D matplotlib

To view the 3D plot in an interactive setting, it is highly recommended that you use Spyder from Anaconda.
In order to generate a 3D plot that you can intereact with, make sure you have the following setting:
Preferences > IPython console > Graphics > Backend (or Graphics Backend) > 'Automatic'
You may need to restart / open a new Ipython console after applying the changes
'''

#-------------- Setting the Size --------------------------------
N_max = 8     
N_min = 2     #Define the min / max graph sizes (not recommended going over N=2 / M=10)
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


#---------------- 3D Plot -----------------------
x = zeros( len(N_axis) * len(M_axis) )                 #creating arrays for 3D points
y = zeros( len(N_axis) * len(M_axis) )
z = zeros( len(N_axis) * len(M_axis) )

for i in arange( len(N_axis) ):
    for j in arange( len(M_axis)):
        x[i*len(N_axis)+j] = N_axis[i]
        y[i*len(N_axis)+j] = M_axis[j]
        z[i*len(N_axis)+j] = Angles[i,j]               #setting z-values as eigenangles

fig = plt.figure(facecolor='white',figsize=(10,6))     #plotting the results
ax = fig.gca(projection='3d')
ax.scatter(x,y,z, linewidth=0.2,color="blue",antialiased=True)
ax.set_zlim(0, max(z))
ax.set_xlim(0, N_max)
ax.set_ylim(0, M_max)
ax.set_xlabel('N')
ax.set_ylabel('M')
ax.set_zlabel('Eigen Angle')
plt.draw()
