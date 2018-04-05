import Ntreefunctions as nt
import matplotlib.pyplot as plt

'''
This code is designed to calculate the moments when the quantum systems have a peak path /peak F probability

The user specifies a specific graph size N and M.

A reccomended number of steps for a good-looking plot is provided by the function Ntree_Path_Optimal(N,M). 
Or, the user can specify any desired number of unitary steps.

The values for the peak probabilities, as well as the number of steps, are printed.  

In addition, a plot for Path Probability vs. Unitary Steps is printed.
'''

#-------------- Setting the Size / Steps --------------------------

N = 7                 #define the size of the graph
M = 8                 #for resaonable runtimes, I do not reccomened going over M=20

steps = round(2 * nt.Ntree_Path_Optimal(N,M) )
#steps = 1234          #For a desired number of steps, remove the # from this line and enter the number of steps


#-------------- Running The Code --------------------------------

Path_vec  = []
F_vec     = []         #creates empty arrays that will hold the values for probabilities / steps
Steps_vec = []

r,t = nt.Unitary_Constants(N,M)
total_states = nt.Total_States(N,M)
SMAT = nt.Ntree_Initialize(total_states,M)           #create the matrix representing the initial state of the system 
nt.Ntree_Unitary_Check(SMAT,N,M)

pp = 0
pf = 0
Max_P = 0
Max_S = 0                                            #create variables to hold information about peak values / steps
Max_Pf = 0
Max_Sf = 0
step = 0


for i in arange(steps):                              #run the desired number of unitary steps on the system
    step = int( i + 1 )
    SMAT = nt.Ntree_Unitary_Step(SMAT,M,N,r,t)
    nt.Ntree_Unitary_Check(SMAT,N,M)
    pp = nt.Ntree_Probability_Path(SMAT,M)           #calculate the probabilities for the path / F
    pf = nt.Ntree_Probability_F(SMAT,M)
    Path_vec.append(pp*.01)                          #add the values for the probabilities into vectors, for plotting later
    F_vec.append(pf*.01)
    Steps_vec.append(step)
    if( pp > Max_P):
        Max_P = pp
        Max_S = step
    if( pf > Max_Pf):                                #store the max probabilities and when they occur
        Max_Pf = pf
        pp_check = pp
        Max_Sf = step
        
        
        
#------------ Displaying The Results---------------------        
        
print("Maximum Path Probability: ",round(Max_P,2),"%, after ",Max_S," Unitary Steps")      
print("Maximum F Probability:     ",round(Max_Pf,2),"%, after ",Max_Sf," Unitary Steps")   
        
      
plt.figure(facecolor='white',figsize=(10,6)) 
plt.plot(Steps_vec,Path_vec,linewidth=2.2,color='black')
plt.plot(Steps_vec,F_vec,linewidth=2.2,color='red')
plt.title( 'Probability of Measuring a Path State',fontsize=20 )
plt.legend(['Correct Path Probability','Correct Final Node Probability'],loc='upper right',prop={'size':14})
plt.axis([0,len(Steps_vec),0,1])
plt.xlabel(' Steps ',fontsize=18)
plt.tick_params(axis='x',labelsize=18,pad = 8)
plt.tick_params(axis='y',labelsize=18,pad = 8)
plt.show()



