import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

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
            irred = list(f[keys[9]])
        Lambdas = list(f[keys[10]])
        selfenergy = list(f[keys[11]])
        
    if only_SE:
        return selfenergy, Lambdas
    else:
        return K1a, K1p, K1t, \
               K2a, K2p, K2t, \
               K3a, K3p, K3t, \
               irred, selfenergy, Lambdas

def rearrange_SE(selfenergy):
    """
    bring selfenergy into appropriate format: 2 x 2 x nw numpy array, where:
    1st index    : 0 = Re(Sigma), 1 = Im(Sigma)
    2nd index iK : 0 = Sigma^R,   1 = Sigma^K
    3rd index iw : frequency index
    """
    nw = len(selfenergy)/2
    SE = np.zeros((2, 2, nw))
    for iK in range(2):
        for iw in range(nw):
            SE[0, iK, iw] = selfenergy[iK*nw+iw][0]
            SE[1, iK, iw] = selfenergy[iK*nw+iw][1]
    return SE

try:
    filename = sys.argv[1]
except IndexError:
    filename = "testfile.h5"
else:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


""" PARAMETERS """

wmin = -20
wmax = 20
re = 0   # 0 = Re(Sigma), 1 = Im(Sigma)
iK = 0   # Keldysh component: 0 = Sigma^R, 1 = Sigma^K

fs = 14  # font size

""" """ 


selfenergy, Lambdas = load_hdf5(filename, True)

nw = len(selfenergy[0])/2
w = np.linspace(wmin, wmax, nw)

lines = []
leg = []
for i in range(7):
    SE = rearrange_SE(selfenergy[i])
    line, = plt.plot(w, SE[re, iK])
    lines.append(line)
    leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
    
plt.xlabel(r'$\omega$', fontsize=fs)
if not re:
    if not iK:
        plt.ylabel(r'$\mathrm{Re}(\Sigma^R(\omega))$', fontsize=fs)
    else:
        plt.ylabel(r'$\mathrm{Re}(\Sigma^K(\omega))$', fontsize=fs)
else:
    if not iK:
        plt.ylabel(r'$\text{Im}(\Sigma^R(\omega))$', fontsize=fs)
    else:
        plt.ylabel(r'$\text{Im}(\Sigma^K(\omega))$', fontsize=fs)            

plt.xticks(fontsize=fs)        
plt.yticks(fontsize=fs)

plt.legend(lines, leg, fontsize=fs, loc='upper right')

plt.show()
    