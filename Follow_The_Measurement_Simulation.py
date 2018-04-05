import Ntreefunctions as nt
import matplotlib.pyplot as plt
import Quantum_Walk_Simulation as qw
from numpy import *
from scipy import *

'''
This code is meant to simulate measurements and calculate average solving speeds for the 'Follow The Measurement' algorithm
This code is only designed to work for N=2

The user spcifies the upper limit M_max as well as the number of simulations per M size

The code calcuates average solving speeds and plots the results
'''
 
def Determine_Smat(M,all_smat):
    '''
    Takes the values from ALL_SMAT and converts them into a single SMAT
    '''
    Smat = zeros(shape=(2*M,3*M))
    Smat = all_smat[int(M-1),0:int(2*M),0:int(3*M)]
    return Smat

def Measurement(n,M,Smat):
    '''
    Simulates a measurement on the quantum system represented by SMAT
    returns an interger 'result' which encodes what kind of state was measured
    returns an integer 'Measued_State' which encodes the layer M' that was measured
    '''
    measure = rand()
    S_value = Smat[2*(M-1),0]**2         + Smat[2*(M-1)+1,0]**2
    F_value = Smat[2*(M-1),3*(M-1)+2]**2 + Smat[2*(M-1)+1,3*(M-1)+2]**2
    Path_value = 0
    for i in arange(M):
        Q = 2*(M-1)+1
        P = 3*i
        if( i == 0):
            Path_value = Path_value + Smat[Q,P+2]**2
        if( i == (M-1) ):
            Path_value = Path_value + Smat[Q,P]**2
        if( 0 < i < (M-1) ):
            Path_value = Path_value + Smat[Q,P]**2 + Smat[Q,P+2]**2
    f_value = 0
    for j in arange(M):
        Q = 2*j
        P = 3*j
        N = ( n**(M - (j+1)) )*(n-1)
        f_value = f_value + N*(Smat[Q,P+1]**2 + Smat[Q+1,P+1]**2)
    M_vec = array([S_value,F_value,Path_value,f_value])
    Measure = 0
    result = 5
    done = False
    for k in arange(4):
        Measure = Measure + M_vec[k]
        if( (Measure >= measure) & (done == False) ):
            result = k+1
            done = True
    if( result == 1):
        Measured_State = M
    if( result == 2):
        Measured_State = 100
    if( result == 3):
        Measure = 0
        path_measure = measure - S_value - F_value
        found_path = False
        Q = int(2.0*(M-1))
        for i in arange(M-1):
            P = 3.0*(i)
            Measure = Measure + Smat[int(Q+1),int(P+2)]**2 + Smat[int(Q+1),int(P+3)]**2
            if( (Measure >= path_measure) & (found_path == False) ):
                found_path = True
                Measured_State = (M-1) - i
    if( result == 4):
        Measured_State = 101
    if( result == 5):
        off_measure = measure - sum(M_vec)
        Measured_State = Off_State_Measure(n,M,off_measure,Smat)
    return result,Measured_State

def Off_State_Measure(n,M,measure,Smat):
    '''
    Used for the function Measurement
    Simulates a measurement on states not along the correct path in the quantum system represented by SMAT
    '''
    Measure = 0
    found = False
    Measured_State = 55
    for j in arange(M-1):
        for i in arange(M-j-1):  
            Q = 2*(j+i+1)
            P = 3*(i)
            N = ( n**(M - (j+2)) )*(n-1)
            Measure = Measure + N*(Smat[Q,P+1]**2 + Smat[Q+1,P+1]**2)
            if( (Measure >= measure) & (found == False) ):
                found = True
                Measured_State = 1 + j
    return Measured_State

def Measurement1(M,Smat):
    '''
    Simulates a measurement on the quantum system represented by SMAT, where M=1
    returns an interger 'result' which encodes what kind of state was measured
    '''
    measure = rand()
    S_value = Smat[0,0]**2 + Smat[1,0]**2
    F_value = Smat[0,2]**2 + Smat[1,2]**2
    f_value = Smat[0,1]**2 + Smat[1,1]**2
    M_vec = array([S_value,F_value,f_value])
    Measure = 0
    done = False
    for k in arange(3):
        Measure = Measure + M_vec[k]
        if( (Measure >= measure) & (done == False) ):
            result = k+1
            done = True
    if( result == 3):
        result = 4
    return result

def Dead_Measurement(n,M):
    '''
    Simulates a measurement on the quantum system represented by a Dead Tree
    returns an interger 'result' which encodes what kind of state was measured
    returns an integer 'Measued_State' which encodes the layer M' that was measured
    '''
    measure = rand()
    Measured_State = 200
    total_states = 0.0
    for j in arange(0,M+1):
        jj = int(j)
        total_states = total_states + n**jj
    S_value = 1.0/total_states
    F_value = (n**M) / total_states
    M_vec = array([S_value,F_value])
    Measure = 0
    result = 3
    done = False
    for k in arange(2):
        Measure = Measure + M_vec[k]
        if( (Measure >= measure) & (done == False) ):
            result = k+1
            done = True
    if( result == 1):
        Measured_State = 201
    if( result == 2):
        Measured_State = 202
    if( result == 3):
        Measure2 = 0
        Miss_Measure = Measure - S_value - F_value
        done2 = False
        for i in arange(M-1):
            Measure2 = Measure2 + (n**(M-1-i))/total_states
            if( (Measure2 >= Miss_Measure) & (done2 == False) ):
                Measured_State = M-1-i
                done2 = True      
    return result,Measured_State

#-------------- Setting the Size --------------------------------------
    
N = 2           #do not change from N=2

M_max = 14       #Define the max M value that will be simulated    (max allowed M is 28)
                 #For reasonable runtimes, I do not recommend exceding M=15

#-------------- Generate all the SMATs --------------------------------

ALL_SMAT = qw.ALL_SMAT_Quantum_Walk(N,M_max)
                                    
                                    
#--------------- Setup Arrays -----------------------                                    
Optimal_Steps = zeros(M_max)        #array for holding optimal number of steps per M size
M_vec         = zeros(M_max)        #arrays for plotting
Limit_vec     = zeros(M_max)
F_vec         = zeros(M_max)
Classical_vec = zeros(M_max)

for i in arange(1,M_max+1):         #set the values in the arrays
    M = int(i)
    M_vec[i-1] = M
    Limit_vec[i-1] = 1.0
    if( M!=1 ):
        Optimal_Steps[i-1] = nt.Ntree_Path_Optimal(N,M)
        F_vec[i-1] = (( (1.0) / nt.Ntree_F_Prob(N,M) ) * nt.Ntree_F_Optimal(N,M) ) / Optimal_Steps[i-1]
    else:
        Optimal_Steps[i-1] = 4
        F_vec[i-1] = 1.0 / ((50.0)/54)
    Classical_vec[i-1] = ( (0.5)*nt.Total_States(N,M) ) / Optimal_Steps[i-1]
  

#------------- Run the Simulations -----------------------------------

Trials = int( 1500 )                          #choose how many simulations to run per size M (reccomended do not exceed 10**4)
Highest_Layer = M_max                         #sets the largest graph size

Speed_Vec = zeros(M_max)

for kk in arange(1,Highest_Layer+1):          #run the simluations 
    Starting_Layer = int(kk)
    Steps_Vec = zeros(Trials)
    for i in arange(Trials):                                              #run the simluation the desired number of trials
        found = False
        layer = Starting_Layer
        measured_layer = 99
        while( found == False):
            Smat = Determine_Smat(layer,ALL_SMAT)
            if( layer != 1):
                result,next_layer = Measurement(N,layer,Smat)            #perform a simulated measurement
            else:
                result = Measurement1(layer,Smat)
            Steps_Vec[i] = Steps_Vec[i] + Optimal_Steps[layer-1]
            if( result == 1):                                            #measurement result is S
                layer = layer
            if( result == 2):                                            #measurement result is F
                found = True
            if( result == 3):                                            #measurement result is on the correct path
                layer = next_layer
            if( result == 4):                                            #measurement result is a wrong final node
                layer = Starting_Layer
            if( result == 5):                                            #measurement result is into a Dead Tree
                dead_bool = True
                layer = next_layer
                while( dead_bool == True ):
                    d_result,d_nextlayer = Dead_Measurement(N,layer)
                    Steps_Vec[i] = Steps_Vec[i] + Optimal_Steps[layer-1]
                    if( d_result == 1):                                   #measurement result is S
                        layer = layer
                    if( d_result == 2):                                  #measurement result is an end node
                        dead_bool = False
                        layer = Starting_Layer                       
                    if( d_result == 3):                                  #measurement result is not an end node
                        if( layer != 1 ):
                            layer = d_nextlayer
                        else:
                            dead_bool = False
                            layer = Starting_Layer
                                         

    Speed_Vec[kk-1] = round((sum(Steps_Vec)/len(Steps_Vec))/Optimal_Steps[Starting_Layer-1],6)      #calculates average solving speeds
    print ("M:",kk,"| Average Solving Speed: ",Speed_Vec[kk-1])                                     #print the results


#----------- Plotting The Results -------------------


plt.figure(facecolor='white',figsize=(10,6)) 
plt.scatter(M_vec,Speed_Vec,color='red',marker='o')
plt.plot(M_vec,Speed_Vec,linewidth=1.0,linestyle='--',color='red')
plt.scatter(M_vec,F_vec,color='green',marker='^')
plt.plot(M_vec,F_vec,linewidth=1.0,linestyle=':',color='green')
plt.scatter(M_vec,Classical_vec,color='black',marker='s')
plt.plot(M_vec,Classical_vec,linewidth=1.0,color='black')
plt.plot(M_vec,Limit_vec,linewidth=.6,linestyle='--',color='blue')
plt.title( 'Average Solving Speeds',fontsize=24 )
plt.legend(['Follow The Measurement','Search for F Directly','Classical Depth-First'],loc='upper left',prop={'size':16})
plt.axis([.9,M_vec[-1]+.1,.5,max(F_vec)*(1.3)])
plt.xlabel(' M ',fontsize=18)
plt.text(round(M_vec[-1]/2.0), 1.05,'Quantum Theoretical Limit',color='blue', fontsize=12)
plt.tick_params(axis='x',labelsize=15,pad = 8)
plt.tick_params(axis='y',labelsize=15,pad = 8)
plt.show() 


