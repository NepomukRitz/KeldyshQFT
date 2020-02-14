import numpy as np

def set_parameters(selfenergy):
    
    """ PARAMETERS """
    
    wmin = -20
    wmax = 20
    
    epsilon = 0
    Gamma = 1
    
    nLambda = 14
    
    fs = 14  # font size
    
    """ """ 
    
    nw = len(selfenergy[0])/2
    w = np.linspace(wmin, wmax, nw)
    
    return w, epsilon, Gamma, nLambda, fs



