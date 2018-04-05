from numpy import *
from scipy import *


'''
This code stores all of the functions called upon by other codes, for handling quantum systems
'''


def Total_States(N,M):
    '''
    Calculates the total number of states for the Hilbert Space
    '''
    total_states = long(0)
    for j in arange(0,M+1):
        jj = long(j)
        total_states = long(total_states + long(N**jj))
    total_states = long(total_states * 2.0)
    return total_states

def Unitary_Constants(N,M):
    '''
    calculates the reflection / transmissions coefficients t and r
    '''
    r = (N-1.0)/(N+1.0)    #Note that these are different from the theoretical values
    t = (2.0)/(N+1.0)      #These take into account that an N-tree has (N+1) connections per node
    return r,t

def Ntree_Initialize(total,M):
    '''
    Creates the matrix SMAT, that will hold the information for all the states in the system
    The SMAT that is returned represents all the states in the system in an equal superposition
    '''
    Mat   = zeros(shape=(2*M,3*M))
    Norm = 1/sqrt(total)
    for i in arange(0,M):
        for j in arange(0,M):
            Q = 2*i
            P = 3*j
            if( (j <= i) & (i != (M-1)) ):
                Mat[Q,P+1]     = Norm
                Mat[Q+1,P+1]   = Norm
            if( i == (M-1)):
                Mat[Q,P+1]     = Norm
                Mat[Q+1,P]     = Norm
                Mat[Q+1,P+1]   = Norm
                Mat[Q+1,P+2]   = Norm
                if( j == 0):
                    Mat[Q,P]   = Norm
                if( j ==(M-1)):
                    Mat[Q,P+2] = Norm
    return Mat

def Ntree_Unitary_Check(Mat,n,M):
    '''
    Sums the probabilities of all the ststes in the system, to make sure that they sum to 1
    If they do not, prints a messege
    '''
    Prob = 0
    for i in arange(0,M):
        ii = long(i)
        N = long(( n**(M - (ii+1)) )*(n-1))
        for j in arange(0,M):
            Q = 2*i
            P = 3*j
            if( (j <= i) & (i != (M-1)) ):
                Prob = Prob + N*( Mat[Q,P+1]**2 + Mat[Q+1,P+1]**2 )
            if( i == (M-1) ):
                Prob = Prob + Mat[Q+1,P]**2 + Mat[Q+1,P+2]**2
                Prob = Prob + N*( Mat[Q,P+1]**2 + Mat[Q+1,P+1]**2 )
                if( j == 0):
                    Prob = Prob + Mat[Q,P]**2
                if( j == (M-1)):
                    Prob = Prob + Mat[Q,P+2]**2
    Prob = round(Prob,8)
    if( Prob != 1.0 ):
        print("NOT UNITARY. Total Probability in the System: ",Prob)

def Ntree_Unitary_Step(Mat,M,n,r,t):
    '''
    Advances the quantum system 1 time-step according to the unitary operations for the scattering quantum random walk
    '''
    Mf = zeros(shape=(2*M,3*M))       #creates an empty matrix to hold all the final values
    Mc = zeros(shape=(2*M,3*M))       #creates an empty matrix that will copy all the values from Mat
    Mc = Mat[:,:]
    for i in arange(0,M):
        for j in arange(0,M):
            Q = 2*i
            P = 3*j
            if( (j < i ) & (i < (M-1) ) ):
                Mf[Q,P+1] = n*t*Mc[Q-2,P+1] - r*Mc[Q+1,P+1]
                Mf[Q+1,P+1] = t*Mc[Q+3,P+1] + (t*(n-1.0)-r)*Mc[Q,P+1]
            if( (i == j) & ( j!=(M-1)) ):
                Mf[Q,P+1] = Mc[Q+1,P+1]
                Mf[Q+1,P+1] = t*Mc[Q+3,P+1] + (t*(n-1.0)-r)*Mc[Q,P+1]
            if( (i==(M-1)) & (j!=0) & (j!=(M-1) )):
                Mf[Q,P+1] =    n*t*Mc[Q-2,P+1]        - r*Mc[Q+1,P+1]
                Mf[Q+1,P+1] =  t*( Mc[Q+1,P-1] + Mc[Q+1,P+3] ) + (t*(n-2.0)-r)*Mc[Q,P+1]
                Mf[Q+1,P] =    (n-1)*t*Mc[Q,P+1] + t*Mc[Q+1,P+3] - r*Mc[Q+1,P-1]
                Mf[Q+1,P+2] =  (n-1)*t*Mc[Q,P+1] + t*Mc[Q+1,P-1] - r*Mc[Q+1,P+3]
            if( (i==(M-1)) & (j==0) ):
                Mf[Q,P] =    (1)*Mc[Q+1,P]     #S flips -- can change to 1 or -1
                Mf[Q,P+1] =    n*t*Mc[Q-2,P+1]        - r*Mc[Q+1,P+1]
                Mf[Q+1,P+1] =  t*( Mc[Q,P] + Mc[Q+1,P+3] ) + (t*(n-2.0)-r)*Mc[Q,P+1]
                Mf[Q+1,P] =    (n-1)*t*Mc[Q,P+1] + t*Mc[Q+1,P+3] - r*Mc[Q,P]
                Mf[Q+1,P+2] =  (n-1)*t*Mc[Q,P+1] + t*Mc[Q,P]     - r*Mc[Q+1,P+3]
            if( (i==(M-1)) & (j==(M-1)) ):
                Mf[Q,P+2] =  (-1)*Mc[Q+1,P+2]   #F flips -- can change to 1 or -1
                Mf[Q,P+1] =   Mc[Q+1,P+1]
                Mf[Q+1,P+1] =  t*( Mc[Q+1,P-1] + Mc[Q,P+2] ) + (t*(n-2.0)-r)*Mc[Q,P+1]
                Mf[Q+1,P] =    (n-1)*t*Mc[Q,P+1] + t*Mc[Q,P+2] - r*Mc[Q+1,P-1]
                Mf[Q+1,P+2] =  (n-1)*t*Mc[Q,P+1] + t*Mc[Q+1,P-1] - r*Mc[Q,P+2]
    return Mf               

def Ntree_Probability_Path(Mat,M):
    '''
    Calculates the total probability accumulated along the states that make up the path from S to F
    '''
    Prob = 0
    for j in arange(1,M):
        Q = 2*(M-1)
        P = 3*j
        Prob = Prob + Mat[Q+1,P]**2 + Mat[Q+1,P+2]**2
        if( j == 1):
            Prob = Prob + Mat[Q+1,P-1]**2        
        if( j == (M-1)):
            Prob = Prob + Mat[Q,P+2]**2
    Prob = round(Prob,7)
    Prob = Prob * 100.0
    return Prob
    
def Ntree_Probability_F(Mat,M):
    '''
    Calculates the total probability accumulated along the states connected to the node F
    '''
    Prob = 0
    Q = 2*(M-1)
    P = 3*(M-1)
    Prob = Prob + Mat[Q+1,P+2]**2 + Mat[Q,P+2]**2
    Prob = round(Prob,7)
    Prob = Prob * 100.0
    return Prob
    
def Ntree_Probability_S(Mat,M):
    '''
    Calculates the total probability accumulated along the states connected to the node S
    '''
    Prob = 0
    Q = 2*(M-1)
    P = 0
    Prob = Prob + Mat[Q,P]**2 + Mat[Q+1,P]**2
    Prob = round(Prob,7)
    Prob = Prob * 100.0
    return Prob

def Ntree_Path_Distribution(M,Smat,Pmat,s):
    '''
    Returns a matrix that contains the probabilities for individual path states
    '''
    step = int(s)
    for j in arange(1,M):
        Q = 2*(M-1)
        P = 3*j
        Pmat[j-1,step] = Smat[Q+1,P]**2 + Smat[Q+1,P-1]**2       
        if( j == (M-1)):
            Pmat[j,step] = Smat[Q,P+2]**2 + Smat[Q+1,P+2]**2
    return Pmat


def Classical_Linear_Approx(N,M):
    '''
    Generates arrays X and Y that represent the linear approximation for the classical search
    '''
    once = True
    total_nodes = int( (0.5)*Total_States(N,M) + 1)
    X = zeros( total_nodes + 1 )
    Y = zeros( total_nodes + 1 )
    for i in arange( total_nodes + 1 ):
        X[i] = i
        Y[i] = float(i) / (total_nodes)
    return X,Y

def Quantum_F_Search(N,M,prob,steps):
    '''
    Generates arrays X and Y that represent probability for the Searching for F Directly search algorithm
    '''    
    total = Total_States(N,M)
    X = zeros( total + 1 )
    Y = zeros( total + 1 )
    step = 0
    F_steps = Ntree_F_Optimal(N,M)
    P_fail = 1.0 - Ntree_F_Prob(N,M)
    pf = 1.0
    for i in arange( total + 1 ):
        X[i] = i
        if( step >= F_steps ):
            pf = pf*P_fail
            step = 0
        Y[i] = 1.0 - pf
        step = step + 1
    return X,Y

#------------------------------------------------------------
#-------              Stored Values Found Numerically
#------------------------------------------------------------
    
#The following functions are all designed to return values for Probabilities / Peak Locations
#All of the following values are exact peaks, found numerically by the code Calculating_Peaks.py
#The main reason for storing the values is to save on computation time in other codes


def Ntree_Path_Optimal(N,M):
    '''
    Returns the optimal number of unitary steps to run a quantum walk, based on N and M, for a peak path measurement
    N = {2,15}   |   M = {2,28}      allowed values
    If an N and M combination is entered, but the value is unknown, returns a 0 or 'out of bounds' error
    The values that are returned were all found numerically
    '''
    Steps_Mat = zeros(shape=(14,27))      #first index = N, second index = M
    Steps_Mat[0,:]    = array([4,11,18,27,42,61,97,145,214,321,476,711,1042,1531,2247,3276,4755,6875,9965,14409,20869,30186,	43584,	63104,	90930,	131056,188933])
    Steps_Mat[1,0:20] = array([8,19,36,65,126,236,439,810,1474,2663,4787,8560,15348,27417,48872,87616,156938,280612,499326,888835])
    Steps_Mat[2,0:17] = array([10,25,59,131,283,603,1282,2705,5721,12023,25277,52597,109403,226895,471112,973240,2007350])
    Steps_Mat[3,0:11] = array([13,35,87,217,528,1263,2998,7201,17093,40267,94550])
    Steps_Mat[4,0:10] = array([14,43,118,323,882,2357,6241,16191,42002,108302])
    Steps_Mat[5,0:9]  = array([16,53,160,481,1398,4053,11420,32315,90109])
    Steps_Mat[6,0:9]  = array([18,63,205,673,2086,6377,19476,58680,173546])
    Steps_Mat[7,0:9]  = array([20,73,266,891,2947,9698,30937,97917,308068])
    Steps_Mat[8,0:9]  = array([24,91,328,1167,4064,13879,46611,154809,517530])
    Steps_Mat[9,0:9]  = array([26,103,390,1479,5403,19133,67341,236910,834410])
    Steps_Mat[10,0:9] = array([28,115,469,1841,6930,25666,95173,352193,1287464])
    Steps_Mat[11,0:9] = array([29,129,554,2245,8659,33726,131326,502612,1911799])
    Steps_Mat[12,0:9] = array([32,147,642,2666,10808,44017,17636,	700011,2773488])
    Steps_Mat[13,0:9] = array([34,163,741,3127,13255,55934,231764,957253,3915556])
    steps = Steps_Mat[N-2,M-2]
    return steps


def Ntree_Path_Prob(N,M):
    '''
    Returns the probability corresponding to the optimal number of unitary steps to run a quantum walk, based on N and M, for a peak path measurement
    N = {2,15}   |   M = {2,28}      allowed values
    If an N and M combination is entered, but the value is unknown, returns a 0 or 'out of bounds' error
    The values that are returned were all found numerically
    '''
    Prob_Mat = zeros(shape=(14,27))      #first index = N, second index = M
    Prob_Mat[0,:]    = array([70.27239,77.28864,73.86883,75.54255,79.55192,82.23506,84.53016,86.79631,88.63398,90.06396,91.27692,92.32139,93.18329,93.93183,94.59822,95.18545,95.66235,96.02988,96.30762,96.51874,96.69035,96.8426	,96.96786,97.09601,97.20999,97.3068,97.39608])
    Prob_Mat[1,0:20] = array([78.92691,81.65839,84.20266,87.82814,89.79587,91.98831,93.773,95.15291,96.16789,96.8308,97.244,97.45043,97.54566,97.5481,97.48112,97.41911,97.42408,97.47998,97.55794,97.64588])
    Prob_Mat[2,0:17] = array([89.11384,86.57506,91.24128,93.99362,95.53879,96.38596,96.77754,96.91084,96.99907,97.06291,97.23608,97.37224,97.49586,97.62794,97.78411,97.94902,98.1151])
    Prob_Mat[3,0:11] = array([93.31686,93.27671,94.44446,95.93268,96.63097,96.93756,96.89019,96.89695,97.16434,97.43302,97.67978])
    Prob_Mat[4,0:10] = array([94.92633,94.69624,96.29646,96.44852,96.38108,96.64287,97.07532,97.43046,97.72561,98.05426])
    Prob_Mat[5,0:9]  = array([94.78716,96.0279,95.84256,96.10377,96.59024,97.15061,97.54492,97.97793,98.36547])
    Prob_Mat[6,0:9]  = array([94.03591,95.19431,95.54122,96.42153,97.0693,97.5737,98.07148,98.55917,98.89551])
    Prob_Mat[7,0:9]  = array([93.09059,94.9381,95.96579,96.8523,97.48599,98.10703,98.62339,98.89598,98.95758])
    Prob_Mat[8,0:9]  = array([93.07506,95.05914,96.50435,97.21374,97.90912,98.58647,98.92936,98.90707,98.85793])
    Prob_Mat[9,0:9]  = array([93.98705,95.72333,96.79383,97.58425,98.48529,98.84533,98.87058,98.83237,98.88328])
    Prob_Mat[10,0:9] = array([94.66728,96.10909,97.03059,98.0798,98.75486,98.88207,98.78156,98.85789,98.97845])
    Prob_Mat[11,0:9] = array([94.96224,96.36036,97.44035,98.55122,98.83314,98.74555,98.79835,98.94838,99.06538])
    Prob_Mat[12,0:9] = array([95.14741,96.55664,97.96027,98.74119,98.75862,98.7306,98.90433,99.03604,99.19117])
    Prob_Mat[13,0:9] = array([95.31697,96.89523,98.36296,98.79555,98.65665,98.80768,98.97569,99.12769,99.3439])
    prob = (.01)*Prob_Mat[N-2,M-2]
    return prob


def Ntree_F_Optimal(N,M):
    '''
    Returns the optimal number of unitary steps to run a quantum walk, based on N and M, for a peak F measurement
    N = {2,15}   |   M = {2,28}      allowed values
    If an N and M combination is entered, but the value is unknown, returns a 0 or 'out of bounds' error
    The values that are returned were all found numerically
    '''
    Steps_Mat = zeros(shape=(14,27))      #first index = N, second index = M
    Steps_Mat[0,:]    = array([3,9,17,31,47,61,83,135,215,289,451,649,929,1721,2259,3429,4935,7009,9921,14193,20199,28655,40949,62653,90107,128149,181897])
    Steps_Mat[1,0:20] = array([11,19,27,45,87,247,325,833,1493,2603,5037,8765,15307,26763,46917,81507,146657,258075,449527,942153])
    Steps_Mat[2,0:17] = array([11,31,59,135,289,585,1251,2509,5149,12489,25311,51163,107507,216431,466967,939261,2150731])
    Steps_Mat[3,0:11] = array([13,35,97,223,549,1285,2871,6529,1461,33711,97293])
    Steps_Mat[4,0:10] = array([15,47,117,299,791,2009,5147,15309,45543,113297])
    Steps_Mat[5,0:9]  = array([15,51,137,493,1313,3639,12291,33069,88531])
    Steps_Mat[6,0:9]  = array([17,55,199,603,2125,7451,20599,61645,176869])
    Steps_Mat[7,0:9]  = array([20,71,223,1035,3181,9921,30319,100787,307401])
    Steps_Mat[8,0:9]  = array([19,75,293,1363,4455,14529,46923,149851,479615])
    Steps_Mat[9,0:9]  = array([19,79,321,1123,5477,18799,68883,233283,777145])
    Steps_Mat[10,0:9] = array([21,145,541,2007,7207,25693,91147,319391,1239209])
    Steps_Mat[11,0:9] = array([35,151,593,2247,8617,31197,125779,462801,2077461])
    Steps_Mat[12,0:9] = array([37,173,701,2549,10845,41725,159135,725447,2942801])
    Steps_Mat[13,0:9] = array([39,179,759,3127,12561,49095,257455,1012879,3960792])
    steps = Steps_Mat[N-2,M-2]
    return steps
 
def Ntree_F_Prob(N,M):
    '''
    Returns the probability corresponding to the optimal number of unitary steps to run a quantum walk, based on N and M, for a peak F measurement
    N = {2,15}   |   M = {2,28}      allowed values
    If an N and M combination is entered, but the value is unknown, returns a 0 or 'out of bounds' error
    The values that are returned were all found numerically
    '''
    Prob_Mat = zeros(shape=(14,27))      #first index = N, second index = M
    Prob_Mat[0,:]    = array([39.68254,39.73114,33.44412,24.721,18.32354,17.05765,14.08095,13.7107,11.59615,10.3428,9.19513,8.18942,7.71781,7.38734,6.96936,6.49453,6.0879,5.7558,5.43753,5.15645,4.86812,4.61873,4.37788,4.18241,4.00891,3.85122,3.69559])
    Prob_Mat[1,0:20] = array([50.04988,38.6847,28.15746,20.98355,17.30811,16.84603,11.95402,12.56653,11.39169,10.05256,9.19668,8.43564,7.79482,7.20519,6.65079,6.17019,5.74967,5.3734,5.0375,4.7855])
    Prob_Mat[2,0:17] = array([66.22468,45.24174,33.60415,26.24908,20.51663,16.43662,13.97285,12.11895,10.58353,9.55563,8.66018,7.93145,7.28255,6.7413,6.26856,5.88603,5.61879])
    Prob_Mat[3,0:11] = array([70.89498,46.39835,33.6597,24.39172,19.49772,16.3192,13.74858,11.81932,10.27686,9.01891,8.44828])
    Prob_Mat[4,0:10] = array([67.95099,43.51865,33.39082,23.93143,19.00505,15.23845,12.46022,11.19698,10.08758,9.28726])
    Prob_Mat[5,0:9]  = array([67.39495,44.83509,29.53937,21.73294,17.65254,14.56908,12.89504,11.54116,10.33849])
    Prob_Mat[6,0:9]  = array([65.3178,42.13036,29.36327,21.57953,16.92662,14.64468,13.05227,11.66564,10.49])
    Prob_Mat[7,0:9]  = array([61.81768,40.28396,27.96121,20.78178,17.63513,15.11988,13.14911,11.54911,10.37373])
    Prob_Mat[8,0:9]  = array([60.10614,39.43954,26.31609,21.13471,17.80339,15.329,13.30514,11.64058,10.28439])
    Prob_Mat[9,0:9]  = array([56.78326,37.27224,25.52158,18.53615,18.06886,15.17411,12.97752,11.41841,10.139])
    Prob_Mat[10,0:9] = array([54.95566,36.03263,27.32877,21.7102,17.97606,15.1484,12.98399,11.25927,9.970])
    Prob_Mat[11,0:9] = array([55.96417,37.41838,28.24303,22.02093,17.93681,14.90846,12.67194,11.07184,10.10133])
    Prob_Mat[12,0:9] = array([58.58321,37.78787,27.74604,21.78037,17.58648,14.75639,12.55877,11.19634,10.1568])
    Prob_Mat[13,0:9] = array([59.99579,39.01861,28.51463,21.90043,17.57225,14.4605,12.5805,11.36419,10.26825])
    prob = (.01)*Prob_Mat[N-2,M-2]
    return prob






