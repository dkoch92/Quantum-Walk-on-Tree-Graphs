import Ntreefunctions as nt
import matplotlib.pyplot as plt
import itertools


'''
This code calcultes the average exiting speeds from a Dead Tree for N=2.  More specifically,
this code exactly calculates all the possible paths to measure and end node, and then 
averages the speeds for all these paths (probability * total number of steps).

The user specfifes the upper limit to M. Do not change N

The code calculates the average exiting speeds, and then plots the results.
'''


#------------------ Functions ------------------------------

def Paths_Generator(n,M):      
    '''
    For a given N and M, calculates all possible paths leading to an end node, starting from the highest layer
    Outputs
    '''                    
    Path_Mat = zeros(shape=(n**(M-1),M+1))                
    Line_Mat = zeros(shape=(M+1,n**(M-1)))
    for k in arange(n**(M-1)):
        Line_Mat[0,k] = M
    for i in arange(M):
        ii = i*1.0
        if(i == 0):
            Line_Mat[1,0] = 0.0
        else:
            Line_Mat[1,(n**(i-1)):(n**i)] = ii
    for j in arange(2,M):
        for l in arange(M):
            if( l == 0):
                Line_Mat[j,0] = 0.0
                Line_Mat[j,1] = 0.0
            else:
                Line_Mat[j,(n**(l-1)):(n**l)] = Line_Mat[j-1,0:(n**(l-1))]
    for m in arange(n**(M-1)):
        for n in arange(M):
            Path_Mat[m,n] = Line_Mat[n,m]
    return Path_Mat

def Row_N(X,Y,n):
    '''
    A mathematical function needed for the Omega Function
    '''
    total = 0
    terms = len(X)
    X_order = zeros(terms)
    for j in arange(terms):
        X_order[j] = j
    Combinations = asarray(list(itertools.combinations(X_order,n)))
    C = shape(Combinations)
    Permutations = ones(shape=(C[0],terms))
    Permutations = Permutations + Permutations
    for i in arange(C[0]):
        for k in arange(C[1]):
            i2 = int(Combinations[i,k])
            Permutations[i,i2] = 1.0
    for l in arange(C[0]):
        M_subtotal = 1.0
        A_subtotal = 0
        for q in arange(terms):
            M_subtotal = M_subtotal * X[q]**( Permutations[l,q] )
            if( Permutations[l,q] == 2.0):
                A_subtotal = A_subtotal + Y[q]**( Permutations[l,q] - 1 )
        total = total + M_subtotal * A_subtotal
    return total

def Omega(X,Y):
    '''
    A Mathematical function needed for calculating the speeds of paths
    '''
    M = len(X)
    A_subtotal = 0.0
    M_subtotal = 1.0
    for i in arange(M):
        A_subtotal = A_subtotal + (-1.0)**int(i) * Row_N(X,Y,i)
        M_subtotal = M_subtotal * ( X[i] - 1 )**2
    total = A_subtotal / M_subtotal
    return total    
            
def Dead_Tree_Speed(n,M,Path_M,Opt):
    '''
    Takes all the functions / Paths_Mat and calculates the average exiting speeds
    '''
    speed = 0
    for i in arange(n**(M-1)):
        check_len_bool = True
        for ii in arange(M+1):
            if( (check_len_bool == True) & ( Path_M[i,ii] == 0) ):
                check_len_bool = False
                Length = int(ii*1.0)
        X = zeros(Length)
        Y = zeros(Length)
        Z = zeros(Length)
        for k in arange(Length):
            X[k] = (n**(Path_M[i,k] + 1)-1)
            Y[k] = 1.0*(Opt[int(Path_M[i,k]-1)])/Opt[M-1]
            if( k != (Length - 1) ):
                Z[k] = ( n**( Path_M[i,k] - Path_M[i,k+1] ) / ( n**( Path_M[i,k]+1 ) -1 ) )
            else:                
                Z[k] = ( n**( Path_M[i,k] - 0.0 ) / ( n**( Path_M[i,k]+1 ) -1 ) )
        speed = speed + Omega(X,Y)*prod(Z)
    return speed

        



#-------------- Setting the Size / Steps --------------------------
N = 2                      #Define the size of the graph      (do not change)
#-----------
M_max = 11                 #For resaonable runtimes, I do not reccomened going over M=12

Optimal = zeros(M_max+1)                                 #Creates a vector that stores all the optimal steps
for i in arange( 2,M_max+2 ): 
    Optimal[i-1] = nt.Ntree_Path_Optimal(N,i)

Optimal[0] = 4.0                                         #Inserts the case M=1 manualy

                        
#--------------  Running The Code ------------------------------------
        
Speeds_vec = zeros(M_max)                                #Creates vectors to store values for plotting
M_vec = zeros(M_max)

for i in arange(1,M_max+1):                              #Calculates all the exit speeds up to M_Max
    M = int(i*1.0)
    M_vec[i-1] = M
    Path_M = Paths_Generator(N,M)
    Speed = Dead_Tree_Speed(N,M,Path_M,Optimal)
    Speeds_vec[i-1] = Speed
    #print ("M = ",int(M),"Dead Tree___Dead Steps: ",round(Speed*Optimal[i-1])  )      #remove #'s to print values
    #print "Average Speed = ",Speed," Steps"  
    #print ("  ")
    
    
#------------ Displaying The Results---------------------        
     
plt.figure(facecolor='white',figsize=(10,6))              #Plots all the results
plt.scatter(M_vec,Speeds_vec,color='black')
plt.plot(M_vec,Speeds_vec,linewidth=1.2,color='red')
plt.title( 'Average Exit Speeds',fontsize=26 )
plt.axis([.9,len(M_vec)+.1,1.0,max(Speeds_vec)+0.1])
plt.xlabel(' M ',fontsize=18)
plt.ylabel(' Speed ',fontsize=18)
plt.tick_params(axis='x',labelsize=18,pad = 8)
plt.tick_params(axis='y',labelsize=18,pad = 8)
plt.show()    
    
    