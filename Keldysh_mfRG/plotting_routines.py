import numpy as np
import matplotlib.pyplot as plt
import sys
from load_data import load_hdf5, rearrange_SE, compute_spectralfunction
from setup_figure import finalize_figure_selfenergy, finalize_figure_spectralfunction

def plot_selfenergy(filename, cnt, fs, *args):
    iK = int(sys.argv[cnt + 1])   # Keldysh component: 0 = Sigma^R, 1 = Sigma^K
    re = int(sys.argv[cnt + 2])   # 0 = Re(Sigma), 1 = Im(Sigma)
    
    selfenergy, Lambdas, parameters = load_hdf5(filename, True)
    
    lines = []
    leg = []
    try:
        args[0]
    except IndexError:
        for i in range(len(Lambdas)):
            SE, w = rearrange_SE(selfenergy[i], parameters)
            line, = plt.plot(w, SE[iK, re])
            lines.append(line)
            leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))          
    else:
        i = args[0]
        SE, w = rearrange_SE(selfenergy[i], parameters)
        line, = plt.plot(w, SE[iK, re])
        lines.append(line)
        leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
    
    finalize_figure_selfenergy(lines, leg, re, iK, fs)
    
def plot_spectralfunction(filename, fs, *args):
    selfenergy, Lambdas, parameters = load_hdf5(filename, True)
    
    lines = []
    leg = []
    try:
        args[0]
    except IndexError:
        for i in range(len(Lambdas)):
            A, w = compute_spectralfunction(selfenergy[i], parameters)
            line, = plt.plot(w, A)
            lines.append(line)
            leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
    else:
        i = args[0]
        A, w = compute_spectralfunction(selfenergy[i], parameters)
        line, = plt.plot(w, A)
        lines.append(line)
        leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
        
    finalize_figure_spectralfunction(lines, leg, fs)



