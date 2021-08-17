import matplotlib.pyplot as plt
import sys

def initialize_figure():
    try:
        int(sys.argv[1])
    except ValueError:
        filename = sys.argv[1]
        cnt = 2
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    else:
        filename = "cmake-build-debug/testfile.h5"
        cnt = 1
    return filename, cnt

def finalize_figure_selfenergy(lines, leg, re, iK, fs):
    
    plt.xlabel(r'$\omega$', fontsize=fs)
    if not re:
        if not iK:
            plt.ylabel(r'$\mathrm{Re}(\Sigma^R(\omega))$', fontsize=fs)
        else:
            plt.ylabel(r'$\mathrm{Re}(\Sigma^K(\omega))$', fontsize=fs)
    else:
        if not iK:
            plt.ylabel(r'$\mathrm{Im}(\Sigma^R(\omega))$', fontsize=fs)
        else:
            plt.ylabel(r'$\mathrm{Im}(\Sigma^K(\omega))$', fontsize=fs)  
    plt.xticks(fontsize=fs)        
    plt.yticks(fontsize=fs)
    
    plt.legend(lines, leg, fontsize=fs, loc='upper right')
    
    plt.show()
    
def finalize_figure_spectralfunction(lines, leg, fs):
    
    plt.xlabel(r'$\omega$', fontsize=fs)
    plt.ylabel(r'$\mathcal{A}(\omega)$', fontsize=fs)
    plt.xticks(fontsize=fs)        
    plt.yticks(fontsize=fs)
    
    plt.legend(lines, leg, fontsize=fs, loc='upper right')
    
    plt.show()



