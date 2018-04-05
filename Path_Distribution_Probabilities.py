import Ntreefunctions as nt
import matplotlib.pyplot as plt

'''
This code shows the breakdown in probabilities that make up the correct path, using the function Ntree_Path_Distribution()

The user specifies a specific graph size N and M.  Ntree_Optimal(N,M) provides the number of steps for the peak probability

A plot is then printed, showing the probabilities of the individual states
'''

#-------------- Setting the Size / Steps --------------------------

N = 2                  #define the size of the graph
M = 14                 #for resaonable runtimes, I do not reccomened going over M=20

steps = round(2 * nt.Ntree_Path_Optimal(N,M) )
#steps = 1234          #For a desired number of steps, remove the # from this line and enter the number of steps



#-------------- Running The Code --------------------------------

Paths_Mat = zeros(shape=(M,int(steps)))                #creates empty matrix to store values for path states' probabilties
Steps_vec = zeros(int(steps))

r,t = nt.Unitary_Constants(N,M)
total_states = nt.Total_States(N,M)
SMAT = nt.Ntree_Initialize(total_states,M)             #create the matrix representing the initial state of the system 
nt.Ntree_Unitary_Check(SMAT,N,M)



for i in arange(steps):                                #run the desired number of unitary steps on the system
    SMAT = nt.Ntree_Unitary_Step(SMAT,M,N,r,t)
    nt.Ntree_Unitary_Check(SMAT,N,M)
    Paths_Mat = nt.Ntree_Path_Distribution(M,SMAT,Paths_Mat,i)          #calculate the probabilities for the path states
    Steps_vec[int(i)] = i + 1

        
        
        
#------------ Displaying The Results---------------------        
     
plt.figure(facecolor='white',figsize=(10,6)) 
for i in arange(M):                                                    #plots the individual states
    if( i!=(M-1) ):                                                    
        lw = 1.3 - round((int(i)*1.0)/(M-1.0),1)                       #states closer to F are thinner lines
        plt.plot(Steps_vec,Paths_Mat[i,:],linewidth=lw,color='black')
    else:
        plt.plot(Steps_vec,Paths_Mat[i,:],linewidth=0.25,color='red')  #plots the edge connected to F in red
plt.title( 'Probabilities of Individual Path States',fontsize=20 )
plt.axis([0,len(Steps_vec),0,max(Paths_Mat[M-1,:])+.01])
plt.xlabel(' Steps ',fontsize=18)
plt.tick_params(axis='x',labelsize=18,pad = 8)
plt.tick_params(axis='y',labelsize=18,pad = 8)
plt.show()