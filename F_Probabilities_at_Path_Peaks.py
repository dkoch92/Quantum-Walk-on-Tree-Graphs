import Ntreefunctions as nt
import matplotlib.pyplot as plt
from Quantum_Walk_Simulation import Quantum_Walk


'''
This code shows the probabilities of measuring F, when the system is prepapred for a peak path measurement

This code is designed for N=2 graphs, and start at M=2

The user specifies M, the largest graph size that will be calculated

A plot is then printed, showing a comparison of F probabilities for the two ways to prepre the system
'''
#------------------- Setting the Size ------------------------
N = 2                              #Choose the N-tree graph to simulate      (recommended do not go over 3 )
M_max = 15                         #Recommended Example:   N=2   M=15
 

#------------------- Running The Code -----------------------                
M_vec = zeros(M_max-1)             #Choose the highest M to calcultate up to (recommended do not go over 15)
for i in arange(2,M_max+1):        #Create an array to store M values
    M_vec[i-2] = int(i)

F_path = zeros(len(M_vec))         #Create arrays to store the probability values
F_f    = zeros(len(M_vec))

for k in arange(2):
    if(k==0):                                                  #Compute the probability of measuring F
        for j in arange( len(M_vec) ):                         #for when the system is prepared for a peak path measurement
            M = int(M_vec[j])
            Steps_path = nt.Ntree_Path_Optimal(N,M)
            SMAT_path = Quantum_Walk(Steps_path,N,M)
            F_path[j] = .01*(nt.Ntree_Probability_F(SMAT_path,M))
    if(k==1):
        for j in arange( len(M_vec) ):                         #Compute the probability of measuring F
            M = int(M_vec[j])                                  #for when the system is prepared for a peak F measurement
            Steps_f = nt.Ntree_F_Optimal(N,M)
            SMAT_f = Quantum_Walk(Steps_f,N,M)
            F_f[j] = .01*(nt.Ntree_Probability_F(SMAT_f,M))


plt.figure(facecolor='white',figsize=(10,6))                               #Plot the results
plt.scatter(M_vec,F_f,linewidth=2.2,color='blue')
plt.plot(M_vec,F_f,linewidth=.3,color='blue')
plt.scatter(M_vec,F_path,linewidth=2.2,color='green',marker='^')
plt.plot(M_vec,F_path,linestyle='--',linewidth=.3,color='green')
plt.title( 'Probing at Different Peaks',fontsize=25 )
plt.legend(['Maximum F Peak','Maximum Path Peak'],loc='upper left',prop={'size':18})
plt.axis([1,M_vec[-1],0,F_f[0]+.15])
plt.xlabel(' M ',fontsize=18)
plt.tick_params(axis='x',labelsize=15,pad = 8)
plt.tick_params(axis='y',labelsize=15,pad = 8)
plt.show() 



