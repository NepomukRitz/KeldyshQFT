import matplotlib.pyplot as plt
from numpy import pi, imag
import sys
from load_data import load_hdf5, rearrange_SE, compute_spectralfunction
from setup_figure import initialize_figure, finalize_figure_selfenergy, finalize_figure_spectralfunction
from parameters import set_parameters

filename, cnt = initialize_figure()

plot = sys.argv[cnt]

def plot_selfenergy(*args):
    iK = int(sys.argv[cnt + 1])   # Keldysh component: 0 = Sigma^R, 1 = Sigma^K
    re = int(sys.argv[cnt + 2])   # 0 = Re(Sigma), 1 = Im(Sigma)
    
    selfenergy, Lambdas = load_hdf5(filename, True)
    w, epsilon, Gamma, nLambda, fs = set_parameters(selfenergy)
    
    lines = []
    leg = []
    try:
        args[0]
    except IndexError:
        for i in range(nLambda):
            SE = rearrange_SE(selfenergy[i])
            line, = plt.plot(w, SE[iK, re])
            lines.append(line)
            leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))          
    else:
        i = args[0]
        SE = rearrange_SE(selfenergy[i])
        line, = plt.plot(w, SE[iK, re])
        lines.append(line)
        leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
    
    finalize_figure_selfenergy(lines, leg, re, iK, fs)
    
def plot_spectralfunction(*args):
    selfenergy, Lambdas = load_hdf5(filename, True)
    w, epsilon, Gamma, nLambda, fs = set_parameters(selfenergy)
    
    lines = []
    leg = []
    try:
        args[0]
    except IndexError:
        for i in range(nLambda):
            A = compute_spectralfunction(selfenergy[i], w, epsilon, Gamma)
            line, = plt.plot(w, A)
            lines.append(line)
            leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
    else:
        i = args[0]
        A = compute_spectralfunction(selfenergy[i], w, epsilon, Gamma)
        line, = plt.plot(w, A)
        lines.append(line)
        leg.append(r'$\Lambda=\ $' + str(round(Lambdas[i], 2)))
        
    finalize_figure_spectralfunction(lines, leg, fs)

if plot == 'selfenergy':
    plot_selfenergy(13)
elif plot == 'spectralfunction':
    plot_spectralfunction(0)

