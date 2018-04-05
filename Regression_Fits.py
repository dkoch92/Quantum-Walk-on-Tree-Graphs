from numpy import *
from scipy import *
import matplotlib.pyplot as plt

'''
This code stores all of the functions called upon by other codes, for handling fitting data
'''
def Linear_Regression(X,Y):
    '''
    Inputs arrays of data points X and Y:    (x,y)
    Outputs the constants A and B that give best fits according to a linear function
    Also outputs the correlation coefficient r
    '''
    x_sum = 0
    for i in arange(len(X)):
        x_sum = x_sum + X[i]
    x_m = ((1.0)*x_sum)/len(X)
    y_sum = 0
    for i in arange(len(Y)):
        y_sum = y_sum + Y[i]
    y_m = ((1.0)*y_sum)/len(Y)
    x2_sum = 0
    for i in arange(len(X)):
        x2_sum = x2_sum + (X[i] - x_m)**2
    y2_sum = 0
    for i in arange(len(Y)):
        y2_sum = y2_sum + (Y[i] - y_m)**2
    xy_sum = 0
    for i in arange(len(Y)):
        xy_sum = xy_sum + (X[i] - x_m)*(Y[i] - y_m)
    Sxx = ((1.0)*x2_sum)/len(X)    
    Syy = ((1.0)*y2_sum)/len(Y)
    Sxy = ((1.0)*xy_sum)/len(Y)
    B = Sxy/Sxx
    A = y_m-x_m*B
    r = (Sxy)/( (Sxx)**(.5) * (Syy)**(.5) )
    return A,B,r    

def Power_Regression(X,Y):
    '''
    Inputs arrays of data points X and Y:    (x,y)
    Outputs the constants A and B that give best fits according to a power function
    Also outputs the correlation coefficient r
    '''
    x_sum = 0
    for i in arange(len(X)):
        x_sum = x_sum + log(X[i])
    x_m = ((1.0)*x_sum)/len(X)
    y_sum = 0
    for i in arange(len(Y)):
        y_sum = y_sum + log(Y[i])
    y_m = ((1.0)*y_sum)/len(Y)
    x2_sum = 0
    for i in arange(len(X)):
        x2_sum = x2_sum + (log(X[i]) - x_m)**2
    y2_sum = 0
    for i in arange(len(Y)):
        y2_sum = y2_sum + (log(Y[i]) - y_m)**2
    xy_sum = 0
    for i in arange(len(Y)):
        xy_sum = xy_sum + (log(X[i]) - x_m)*(log(Y[i]) - y_m)
    Sxx = ((1.0)*x2_sum)/len(X)    
    Syy = ((1.0)*y2_sum)/len(Y)
    Sxy = ((1.0)*xy_sum)/len(Y)
    B = Sxy/Sxx
    A = e**(y_m-x_m*B)
    r = (Sxy)/( (Sxx)**(.5) * (Syy)**(.5) )
    return A,B,r

def Power_Fit_Plot(M,X,Y,A,B,r,choice):
    '''
    Takes in data for plotting eigenangles
    Outputs a graph that shows these eigenangles along with best fits
    '''
    plt.figure(facecolor='white',figsize=(10,6))
    plt.title( 'Constant N, Varying M',fontsize=25 )
    N = 10000
    for i in arange(len(M)):
        plot_bool = False
        for b in arange( len(choice)):
            if( choice[b] == M[i] ):
                plot_bool = True
        if( plot_bool == True):
            X_fit = zeros(N)
            Y_fit = zeros(N)
            X_range = X[-1] - X[0] + 2
            for k in arange(2,N):
                n = ((int(k)*1.0)/N)*X_range
                X_fit[k] = n
                Y_fit[k] = A[i]*n**(B[i])
            plt.scatter(X,Y[:,i],color="black")
            plt.plot(X_fit,Y_fit,color="red",linewidth=1.2)
            plot_bool = False
    plt.axis([0,X_range,0,20])  
    plt.xlabel(' M ',fontsize=18)
    plt.show()
    
    
    
    
    
    
    
    
