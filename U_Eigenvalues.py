import Ntreefunctions as nt
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy import *

'''
This code is used for generating Matrices that represent N-Tree Graphs
This code stores all of the functions called upon by other codes, for handling Eigenvalues
'''

#--------------------------------- Functions -----------------

def total_states(M):
    t = 0
    for i in arange(M+1):
        t = t + 2**i
    T = 2.0*t
    return T

def total_states_N(M,N):
    t = 0
    for i in arange(M+1):
        t = t + N**i
    T = 2.0*t
    return T

def Unitary_Check(U,Initial):
    Ut = transpose(U)
    UU = dot(U,Ut)
    for i in arange(len(Initial)):
        for j in arange(len(Initial)):
            if( i == j ):
                if ( round(UU[i,j],10) != 1):
                    print( "not unitary",i,j)
            if( i != j):
                if (round(UU[i,j],10) != 0):
                    print( "not unitary",i,j)
    
            

def Unitary_Eigenvalues(v):
    for i in arange(len(v[0])):
        if( round(abs(v[0][i]),10) != 1.0):
            print( "NOT UNITARY EIGENVALUES" )
            
            
def Orthogonality_Check(v):
    for i in arange(len(v[0])):
        for j in arange(len(v[0])):
            if( i != j):
                c =  around(vdot(v[1][:,i],v[1][:,j]),5)
                #if ( abs(c) != 0.0):
                    #print "not orthogonal",i,j
                    
def Eigenvalues_Check(U,v):
    for i in arange(len(v[0])):
        w = dot(U,v[1][:,i]) - v[0][i]*v[1][:,i]
        for k in arange(len(w)):
            ww = around(w[k],7)
            if( ww != 0.0):
                print ("Eigenvector problem" )                 
                    
def Overlaps(Initial,v):
    B = zeros(len(Initial),dtype=complex)
    for i in arange(len(B)):
        B[i] = vdot(v[1][:,i],Initial)
    s = 0
    p = 0
    for i in arange(len(Initial)):
        s = s + Initial[i]**2
        p = p + abs(B[i])**2
    if( round(s,7) != 1.0 ):
        print ("Initial State Not Unitary",round(s,7))
    if( round(p,7) != 1.0):
        print ("Overlaps Don't Sum to 1",round(p,7))
    B_r = zeros(len(B))
    B_angle = zeros(len(B))
    for i in arange( len(B)):
        B_r[i] = abs(B[i])
        B_angle[i] = angle(B[i]*180/pi)
    return B,B_r,B_angle
                    
def Eigenstates_Polar(v):
    v_r = zeros(shape=(len(v[:,0]),len(v[0,:])))
    v_angle = zeros(shape=(len(v[:,0]),len(v[0,:])))
    for i in arange(len(v[0,:])):
        for j in arange(len(v[:,0])):
            v_r[i,j] = abs(v[i,j])
            v_angle[i,j] = angle(v[i,j])*180/pi
    return v_r,v_angle
            
def Eigenvalues_Polar(v):
    v_angle = zeros(len(v))
    for i in arange(len(v)):
        v_angle[i] = angle(v[i])*180/pi
    return v_angle
                                   
def Display_Polars(vr,va,path,u):
    for i in arange(len(u)):
        for j in arange(len(path)):
            print ("u:",u[i]," state:",path[j],"    r:",vr[j,i],"  A:",va[j,i])

def Find_Largest(v):
    V = zeros(len(v))
    for i in arange( len(v) ):
        V[i] = v[i]
    top2 = zeros(2)
    v_sorted = sort(V)
    for i in arange( len(top2) ):
        found = False
        j = 0
        value = v_sorted[-i-1]
        while ( found == False):
            if( V[j] == value ):
                top2[i] = j
                V[j] = 99999
                found = True
            else:
                j = j + 1
    return top2                                                                                                                                                

def Unitary_Step(U,Psi,steps):
    P = zeros(len(Psi))
    for j in arange(len(Psi)):
        P[j] = Psi[j]
    for i in arange(steps):
        P = dot( U,P )
    return P

def Eigen_Step(Initial,v,BB,steps):
    Psi = zeros(len(Initial))
    for i in arange( len(Initial)):
        Amp = 0
        for j in arange( len(v[1][0,:])):
            Amp = Amp + BB[j]*((v[0][j])**steps)*v[1][i,j]
        if( (Amp) != (conjugate(Amp)) ):
            print ("imaginary amplitude!")
        else:
            Psi[i] = Amp.real
    return Psi
    
def Select_Eigen_Step(Initial,v,BB,steps,states):
    Psi = zeros(len(Initial))
    for i in arange( len(Initial)):
        Amp = 0
        for j in arange( len(states) ):
            jj = states[j]
            Amp = Amp + BB[jj]*((v[0][jj])**steps)*v[1][i,jj]
        if( (Amp) != (conjugate(Amp)) ):
            print ("imaginary amplitude!")
        else:
            Psi[i] = Amp.real
    return Psi
        
def Path_Probability(Psi,states):
    Prob = 0
    for i in arange(len(states)):
        Prob = Prob + abs(Psi[ states[i] ])**2     
    return Prob  
    
        
def Set(H,Int2,Int3,Int4,Int5,Int9,v2,v3,v4,v5,v9,u2,u3,u4,u5,u9):
    path = zeros(2*H)
    for i in arange(H):
        if( i == 0 ):
            path[0] = 2
            path[1] = 3
        else:
            path[2*i]   = path[2*i -2] + 2*(i+1)
            path[2*i+1] = path[2*i -1] + 2*(i+1)
    if( H == 2):
        U = u2
        v = v2
        I = Int2        
    if( H == 3):
        U = u3
        v = v3
        I = Int3   
    if( H == 4):
        U = u4
        v = v4
        I = Int4      
    if( H == 5):
        U = u5
        v = v5
        I = Int5   
    if( H == 9):
        U = u9
        v = v9
        I = Int9
    return U,v,I,path                  
 
#-------------------------------------------------
       
def initial_psi(N,M):
    '''
    Creates an array that represents the initial state of the system
    '''
    states2 = nt.Total_States(N,M)
    Psi = []
    pattern = []
    for j in arange(M+1):
        if( j == 0):
            pattern.append(1.0)
            pattern.append(1.0)
        else:
            pattern.append( sqrt(N-1)*(sqrt(N))**(j-1) )
            pattern.append( sqrt(N-1)*(sqrt(N))**(j-1) )
        for k in arange( len(pattern) ):
            Psi.append( pattern[k] )
    for q in arange( len(Psi) ):
        Psi[q] = (1.0/sqrt(long(states2)))*Psi[q]
    return Psi
                
def Unitary_Matrix(N,M,R,T):
    '''
    Creates the unitary matrix U, that advances the quantum system one time-step
    '''
    R = -R
    l = 0
    for i in arange(1,M+2):
        l = l + 2*i
    U = zeros(shape=(l,l))
    U[0,1] = 1
    for i in arange( l/2.0 ):              #-r / 1 values along diagonal
        k = int(i)
        if( i < (l-2*(M+1))/2.0 ):
            U[2*k+1,2*k] = R
        else:
            if( i == ((l-2*(M+1))/2.0) ):
                U[2*k+1,2*k] = -1.0
            else:
                U[2*k+1,2*k] = 1.0
    q = 0
    p = 0 
    for ii in arange(1,M+1): 
        kk = int(ii)              # t / sqrt(N-1) / r values in M intervals
        q = q + 2*kk
        p = p + 2*(kk-1)
        U[q,p] = T
        U[q+2,p] = sqrt(N-1)*T
        U[q+2,p+(2*kk+1)] = sqrt(N-1)*T
        U[p+1,q+1] = T 
        U[p+1+(2*kk-1),q+1] = R    
        U[p+1,q+3] = sqrt(N-1)*T
        U[p+1+(2*kk-1),q+3] = sqrt(N-1)*T
        U[p+1+(2*kk-1)+2,q+3] = (N-2)*T+R
    q = 6
    p = 2
    for iii in arange(2,M+1):             # increasing t / (n-1)T-R values
        kkk = int(iii)
        q = q + 2*kkk
        p = p + 2*(kkk-1)
        x=-2
        y=-2
        for j in arange( iii-1 ):
            x = x + 2
            y = y + 2
            U[q+x,p+y] = sqrt(N)*T
            U[p+1+y,q+1+x] = sqrt(N)*T
            U[p+1+y+(2*kkk+1),q+1+x] = (N-1)*T+R
    return U

def Find_Eigen(N,M):
    '''
    Input the size of the N-tree graph: N,M
    Outputs a numpy.array that stores all the eigenvalues / eigenvectors
    '''
    r,t = nt.Unitary_Constants(N,M)                     
    U = Unitary_Matrix(N,M,r,t)
    Int = initial_psi(N,M)     
    eigen = LA.eig( U )    
                            
    Unitary_Check( U,Int )                              
    Unitary_Eigenvalues(eigen)
    Orthogonality_Check(eigen)
    Eigenvalues_Check(U,eigen)
    return Int,eigen

def Path_States_Z(M,index,eigen_r,eigen_t):
    '''
    Takes in the eigenvectors / eigenvalues from U, and the index for the leading order
    Outputs the leading order values / eigenangles for the path states
    '''
    i = 2
    j = 4
    eigen_r_path = zeros( 2*M )
    eigen_t_path = zeros( 2*M )
    for k in arange(M):
        eigen_r_path[int(2*k)]   = eigen_r[int(i),index]
        eigen_r_path[int(2*k+1)] = eigen_r[int(i+1),index]
        eigen_t_path[int(2*k)]   = eigen_t[int(i),index]   * pi / 180.0
        eigen_t_path[int(2*k+1)] = eigen_t[int(i+1),index] * pi / 180.0
        i = i + j
        j = j + 2
    return eigen_r_path,eigen_t_path
        
def Path_indices(M):
    '''
    returns an array that stores the index locations for all the path states
    '''
    indices = zeros(2*M)
    i = 2
    j = 4
    for k in arange( M ):
        indices[2*k]   = i
        indices[2*k+1] = i+1   
        i = i + j
        j = j + 2
    return indices   
    
def Eigenangle_Function(N,M):
    '''
    returns the eigenangle predicted by epsilon(N,M) using the most accurate versions of alpha,beta,gamme,rho
    '''
    e = ((47.8743628487)*(M)**(-.551332856993))*(N)**(.0766542527 - .4975196235*M)
    return e






