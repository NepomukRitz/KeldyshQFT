import numpy as np
from numpy import pi, imag
import h5py
from parameters import set_parameters

def load_hdf5(filename, only_SE):
    """
    load hdf5 file and write the data into lists
    
    Parameters:
        filename
        only_SE  : boolean, if true only load selfenergy, otherwise also load 
                   vertex components (# of return arguments depends on this)
    """
    with h5py.File(filename, 'r') as f:
        keys = list(f.keys());
        
        if not only_SE:
            K1a = list(f[keys[0]])
            K1p = list(f[keys[1]])
            K1t = list(f[keys[2]])
            K2a = list(f[keys[3]])
            K2p = list(f[keys[4]])
            K2t = list(f[keys[5]])
            K3a = list(f[keys[6]])
            K3p = list(f[keys[7]])
            K3t = list(f[keys[8]])
            irred = list(f[keys[-4]])
        Lambdas = list(f[keys[-3]])
        parameters = list(f[keys[-2]])
        selfenergy = list(f[keys[-1]])
        
    if only_SE:
        return selfenergy, Lambdas, parameters
    else:
        return K1a, K1p, K1t, \
               K2a, K2p, K2t, \
               K3a, K3p, K3t, \
               irred, selfenergy, Lambdas, parameters

def rearrange_SE(selfenergy):
    """
    bring selfenergy into appropriate format: 2 x 2 x nw numpy array, where:
    1st index iK : 0 = Sigma^R,   1 = Sigma^K
    2nd index    : 0 = Re(Sigma), 1 = Im(Sigma)
    3rd index iw : frequency index
    """
    nw = len(selfenergy)/2
    SE = np.zeros((2, 2, nw))
    for iK in range(2):
        for iw in range(nw):
            SE[iK, 0, iw] = selfenergy[iK*nw+iw][0]
            SE[iK, 1, iw] = selfenergy[iK*nw+iw][1]
    return SE

def compute_spectralfunction(selfenergy, w, epsilon, Gamma):
    """
    compute spectral function from selfenergy
    
    Parameters:
        selfenergy : selfenergy as obtained from rearrange_SE
        w          : list of frequencies
        epsilon    : impurity onsite energy
        Gamma      : impurity-bath coupling
    """
    SE = rearrange_SE(selfenergy)
    SE = SE[0, 0] + 1j*SE[0, 1]
    G = 1/(w - epsilon + 1j*Gamma/2 - SE)
    A = -1/pi * imag(G)
    return A
