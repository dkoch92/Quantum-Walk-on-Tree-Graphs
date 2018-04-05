from Ntreefunctions import *
from numpy import *


def Quantum_Walk(steps,N,M):
    '''
    Simulates the quantum random walk for the desirned number of unitary steps
    Returns SMAT, the matrix holding the values for all the states in the system
    This function is designed to be used for cases where one only wants the state of the system after a desired number of steps
    '''
    r,t = Unitary_Constants(N,M)
    total_states = Total_States(N,M)
    sMAT = Ntree_Initialize(total_states,M)
    Ntree_Unitary_Check(sMAT,N,M)
    for i in arange(steps):
        sMAT = Ntree_Unitary_Step(sMAT,M,N,r,t)
        Ntree_Unitary_Check(sMAT,N,M)
    return sMAT

def ALL_SMAT_Quantum_Walk(N,M_max):
    '''
    Simulates the quantum random walk for the number of unitary steps given by Ntree_Path_Optimal(N,M), for numerous sizes M
    Returns All_SMAT, the matrix holding the values for all the states in the system, for each quantum system
    This function is designed to be used for cases where one only wants the state of the system after a desired number of steps
    '''
    all_smat = zeros(shape=(M_max,2*M_max,3*M_max))
    M= 1
    Smat1 = array([[1.0/3,1.0/3,-5.0/3],[-1.0/3,-1.0/3,5.0/3]])       #Hand M=1 case manually
    Smat1 = Smat1 * (1.0/sqrt(6.0))  
    all_smat[int(M-1),0:int(2*M),0:int(3*M)] = Smat1
    for i in arange(2,M_max + 1):
        M = int(i)
        steps = Ntree_Path_Optimal(N,M)
        smat = Quantum_Walk(steps,N,M)
        all_smat[int(M-1),0:int(2*M),0:int(3*M)] = smat
    return all_smat
        
    