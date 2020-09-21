import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from load_data import load_data, load_SIAM_NRG

# use latex and serif fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def generate_filename_fRG(path, diagClass, nLoops, n1=301, sf=False, n2=201, Gamma=0.2, V=0.0, T=0.01, L_ini=20, nODE=50):
    
    Gamma_str = "%.6f" % Gamma
    V_str = "%.6f" % V
    T_str = "%.6f" % T
    
    filename = "K" + str(diagClass) + "_" + str(nLoops) + "LF_n1=" + str(n1)+"_"
    if(sf):
        filename += "static_"
    elif(diagClass>1):
        filename += "n2=" + str(n2) + "_"
        
    filename += "Gamma=" + Gamma_str + "_"
    
    if(V!=0.0):
        filename += "V=" + V_str + "_"
    if(T!= 0.01):
        filename += "T=" + T_str + "_"
        
    filename += "L_ini=" + str(L_ini) + "_" + "nODE=" + str(nODE) + ".h5"
    
    return path + filename

def change_prefix(path, prefix, filename):
    s = list(filename)
    s.insert(len(path), prefix)
    filename = "".join(s)
    return filename


def plot_fRG(typ, iK, filenames, labels, NRG_info):
    U_NRG = NRG_info[0]
    path_NRG = NRG_info[1]
    
    # load NRG reference data
    Gamma_NRG, w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im = load_SIAM_NRG(U_NRG, path_NRG)
    Delta_NRG = Gamma_NRG/2.
    
    # load fRG data
    DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, Lambda, Delta, U, V, \
    v, Sigma, A, \
    w, K1a, K1p, K1t, \
    w2, v2, K2a, K2p, K2t \
        = load_data(U_NRG, Delta_NRG, filenames)
    
    fs = 18 # font size

    #Plot spectral funciton
    if(typ[0]=="A"):        
        fig, ax = plt.subplots(figsize=(18, 6))

        for i in range(len(filenames)):
            ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i])
            
        ax.plot(w_NRG/Delta_NRG, Aimp*np.pi*Delta_NRG, 'k:', label='NRG')    
        ax.set_ylabel(r'$(\mathcal{A}(\nu)\pi\Delta$', fontsize=fs)
        
        ax.set_xlim(-5, 5)
        ax.set_xlabel(r'$\nu/\Delta$', fontsize=fs)
        ax.tick_params(axis='both', labelsize=fs-1)
    
        ax.legend(fontsize=fs)
    
    #Plot self energy, component iK
    elif(typ[0:4] == "self"):
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))

        for i in range(len(filenames)):
            ax[0].plot(v[i]/Delta[i], Sigma[i][iK].real/Delta[i], '.-', label=labels[i])
            ax[1].plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i])

        if(iK==0):
            ax[0].plot(w_NRG / Delta_NRG, (SE_re - U_NRG / 2) / Delta_NRG, 'k:', label='NRG')
            ax[1].plot(w_NRG / Delta_NRG, SE_im / Delta_NRG, 'k:', label='NRG')

            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^R-U/2)/\Delta$', fontsize=fs)
            ax[1].set_ylabel('$(\mathrm{Im}\Sigma^R)/\Delta$', fontsize=fs)
        else:
            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^K)/\Delta$', fontsize=fs)
            ax[1].set_ylabel('$(\mathrm{Im}\Sigma^K)/\Delta$', fontsize=fs)

        for i in range(2):
            ax[i].set_xlim(-20, 20)
            ax[i].set_xlabel(r'$\nu/\Delta$', fontsize=fs)
            ax[i].tick_params(axis='both', labelsize=fs-1)

        ax[0].legend(fontsize=fs)
        fig.tight_layout(pad=3.0)
    
    #Plot susceptibilities
    elif(typ[0:7] == "susc_sp"):    #spin susceptibility
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))

        for i in range(len(filenames)):
            # real part of spin susceptibility (retarded comp. of K1a)
            ax[0].plot(w[i] / Delta[i], K1a[i][0].real / Delta[i], '.-', label=labels[i])
            # imag part of spin susceptibility (retarded comp. of K1a = -advanced comp. of K1a)
            ax[1].plot(w[i] / Delta[i], -K1a[i][0].imag / Delta[i], '.-', label=labels[i])

        ax[0].plot(w_NRG / Delta_NRG, K1a_re * (U_NRG/Delta_NRG)**2 / 2, 'k:', label='NRG')
        ax[1].plot(w_NRG / Delta_NRG, K1a_im * (U_NRG/Delta_NRG)**2 / 2, 'k:', label='NRG')

        ax[0].set_ylabel(r'$\mathrm{Re}\chi_\mathrm{sp}/\Delta$', fontsize=fs)
        ax[1].set_ylabel(r'$\mathrm{Im}\chi_\mathrm{sp}/\Delta$', fontsize=fs)

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs)
            ax[i].set_xlim(-5,5)
            ax[i].tick_params(axis='both', labelsize=fs -1)
        ax[0].legend(fontsize=fs)
        fig.tight_layout(pad=3.0)
        
       
    elif(typ[0:7] == "susc_ch"):    #charge suscpeptibility
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))

        for i in range(len(filenames)):
            # real part of spin susceptibility (retarded comp. of K1a)
            ax[0].plot(w[i] / Delta[i], (2*K1t[i][0]-K1a[i][0]).real / Delta[i], '.-', label=labels[i])
            # imag part of spin susceptibility (retarded comp. of K1a = -advanced comp. of K1a)
            ax[1].plot(w[i] / Delta[i], -(2*K1t[i][0]-K1a[i][0]).imag / Delta[i], '.-', label=labels[i])

        ax[0].plot(w_NRG / Delta_NRG, (2*K1t_re - K1a_re) * (U_NRG/Delta_NRG)**2 / 2, 'k:', label='NRG')
        ax[1].plot(w_NRG / Delta_NRG, (2*K1t_im - K1a_im) * (U_NRG/Delta_NRG)**2 / 2, 'k:', label='NRG')

        ax[0].set_ylabel(r'$\mathrm{Re}\chi_\mathrm{ch}/\Delta$', fontsize=fs)
        ax[1].set_ylabel(r'$\mathrm{Im}\chi_\mathrm{ch}/\Delta$', fontsize=fs)

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs)
            ax[i].set_xlim(-10,10)
            ax[i].tick_params(axis='both', labelsize=fs -1)
        ax[0].legend(fontsize=fs)
        fig.tight_layout(pad=3.0)
            
            
    return fig, ax
