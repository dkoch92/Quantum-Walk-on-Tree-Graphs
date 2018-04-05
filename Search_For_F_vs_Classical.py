import Ntreefunctions as nt
import matplotlib.pyplot as plt
from Quantum_Walk_Simulation import Quantum_Walk

'''
This code is designed to showcase the difference in probabilities for Classical vs. Quantum

The user specifies a specific graph size N and M.

A reccomended number of steps for a good-looking plot is provided by the function Ntree_F_Optimal(N,M). 
An alternative is provided by Ntree_F_Prob(N,M) for a faster runtime

The code generates the probability distributions for the "Search for F Directly" vs Classical algorithms.  

The code prints the average speeds for each algorithm
A plot of these probability distributions is plotted
'''


#-------------- Setting the Size / Steps --------------------------

N = 2                  #define the size of the graph
M = 16                #for resaonable runtimes, I do not reccomened going over N=10 | M=20



#-------------- Running The Code --------------------------------

F_steps = nt.Ntree_F_Optimal(N,M)                        #Calculating the number of steps for an optimal F measurement

#SMAT = Quantum_Walk(F_steps,N,M)                        #These 2 lines of code represent the slow / more precise method  
#F_Prob = nt.Ntree_Probability_F(SMAT,M)

F_Prob = nt.Ntree_F_Prob(N,M)                                     #For faster runtimes / slightly less accurate results

classical_x,classical_y = nt.Classical_Linear_Approx(N,M)         #Generate arrays for plotting probability distributions
quantum_x,quantum_y = nt.Quantum_F_Search(N,M,F_Prob,F_steps)

avg_classical_steps = int(round( (0.5)*nt.Total_States(N,M) ))         #Calculate the average number of steps needed by each algorithm
avg_quantum_steps = int(round( ( 1.0 / nt.Ntree_F_Prob(N,M) ) * nt.Ntree_F_Optimal(N,M) ))
   
#------------ Displaying The Results---------------------        
        
print("Average Classical Speed       : ",avg_classical_steps," steps ")      
print("Average Quantum Speed         : ",avg_quantum_steps," steps ")   
print("Quantum Speedup over Classical: ",round((1.0)*avg_classical_steps/avg_quantum_steps),2)
        
      
plt.figure(facecolor='white',figsize=(10,6)) 
plt.plot(classical_x,classical_y,linewidth=1.5,linestyle='--',color='black')
plt.plot(quantum_x,quantum_y,linewidth=2.2,color='red')
plt.title( 'Comparison of Probability Distributions',fontsize=28 )
plt.legend(['Classical Search Algorithm','Searching for F Directly (quantum)'],loc='upper left',prop={'size':12})
plt.axis([0,len(classical_x),0,1])
plt.xlabel(' Steps ',fontsize=12)
plt.tick_params(axis='x',labelsize=11,pad = 8)
plt.tick_params(axis='y',labelsize=11,pad = 8)
plt.show()



