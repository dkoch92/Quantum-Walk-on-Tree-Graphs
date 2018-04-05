import Ntreefunctions as nt
import matplotlib.pyplot as plt
import U_Eigenvalues as egn
from numpy import *

'''
This code is used for showcasing the accuracy of using the leading order Beta's

The user specifies the graph size: N,M

The code generates the matrix U, finds it's eigenvalues/eigenvectors via U_Eigenvalues functions, 
and plots the results, comparing the leading order approximation to the exact values
'''

#-------------- Setting the Size --------------------------------

N = 2                  #define the size of the graph (not recommended going over N=2
M = 12                 #for resaonable runtimes, I do not reccomened going over M=15
Steps = int(round(2 * nt.Ntree_Path_Optimal(N,M) ))   #calculates teh number of steps for plotting

#-------------- Running The Code --------------------------------

Int,eigen = egn.Find_Eigen(N,M)

B,B_r,B_angle = egn.Overlaps(Int,eigen)                     #stores the overlaps of each eiegenstate with the initial state of the system
Estates_r,Estates_angle = egn.Eigenstates_Polar(eigen[1])           #stores all the eiegenstates as polar values
Evalues_angle           = egn.Eigenvalues_Polar(eigen[0])           #stores all the eigenvalues as polar values

top2 = egn.Find_Largest(B_r)                                #finds the lcoation of the largest two Beta's
index = int(top2[0])

path_indices = egn.Path_indices(M)
Zr_path,Zt_path = egn.Path_States_Z(M,index,Estates_r,Estates_angle)    #creates arrays to hold path states values

Steps_vec   = zeros(Steps)                              #creates arrays to hold values for plotting
Path_Approx = zeros(Steps)
Path_Actual = zeros(Steps)


for i in arange( Steps ):                               #creates the plots
    Steps_vec[int(i)] = i
    prob = 0
    prob2 = 0
    for j in arange( len(Zr_path)):                     #using leading order approximation
        amp = B_r[index] * Zr_path[j] * 2.0 * cos( B_angle[index]*pi/180.0 + Zt_path[j] + Evalues_angle[index]*pi/180.0 * i )
        prob = prob + amp**2
    for k in arange( len(path_indices) ): 
        amp2 = 0 
        for l in arange( len(B_r) ):                    #use exact values
            amp2 = amp2 + B[l]*eigen[1][int(path_indices[k]),l]*(eigen[0][l])**i
        prob2 = prob2 + abs(amp2)**2
    Path_Approx[int(i)] = prob
    Path_Actual[int(i)] = prob2


#------------------ Print Results -----------------------------

plt.figure(facecolor='white',figsize=(10,6)) 
plt.plot(Steps_vec,Path_Actual,linewidth=1.5,linestyle='-',color='black')
plt.plot(Steps_vec,Path_Approx,linewidth=2.2,linestyle='--',color='blue')
plt.title( 'Analytical vs Exact Path Probability',fontsize=28 )
plt.legend(['Full Solution','Approximate Solution'],loc='upper left',prop={'size':16})
plt.axis([0,Steps,0,1.03])
plt.xlabel(' Steps ',fontsize=18)
plt.tick_params(axis='x',labelsize=15,pad = 8)
plt.tick_params(axis='y',labelsize=15,pad = 8)
plt.show() 

