import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from load_data import load_data, load_SIAM_NRG

# use latex and serif fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def effec_distribution(v, T, stat):
    #stat is either =-1 for fermionic statistics or =+1 for bosonic statistics
    res = np.power(np.tanh(v/(2.*T)), -stat)
    if stat == +1:
        try:
            infinite_index = np.where(res == np.inf)[0][0]
            res[infinite_index] = 0
        except IndexError:
            infinite_index = int(len(v)/2)
            res[infinite_index] = 0
    return res

def generate_Keldysh_component(Im_X, v, T, stat):
    #X is retarded component of object
    return 2.*Im_X*effec_distribution(v, T, stat)
    

def generate_filename_fRG(path, diagClass, nLoops, n1=301, sf=False, n2=201, Gamma=1./3., V=0.0, T=0.01, L_ini=1000, nODE=50):
    
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

def generate_label(path_fRG, filenames):
    labels = [] 
    
    for filename in filenames:
        diag_class = filename[len(path_fRG)+1]
        label = r'$\mathcal{K}_' + diag_class
        if (filename[len(path_fRG)+14 : len(path_fRG)+20] == 'static'):
            label += '\mathrm{-sf}'
        
        elif(diag_class == str(2)):
            label += '\mathrm{-' + filename[len(path_fRG)+3] + 'L}'
            
        label += '$'
        labels.append(label)
    return labels


def plot_fRG(typ, iK, filenames, labels, NRG_info, inset):
    U_NRG = NRG_info[0]
    path_NRG = NRG_info[1]
    possible_NRG = NRG_info[2]
    
    U_to_ins_ylim  = {0.1: -0.0005, 0.5:-0.01, 1:-0.003, 1.5: 0.01, 2:-0.03, 2.5:-0.1, 3:0.5}
    
    if possible_NRG:
        # load NRG reference data
        try:
            plot_NRG=True
            Gamma_NRG, w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im = load_SIAM_NRG(U_NRG,path_NRG)
            Delta_NRG = Gamma_NRG/2.
        
        except KeyError:
            plot_NRG=False
            Delta_NRG=0.5
    else:
        plot_NRG=False
        Delta_NRG=0.5
    
    
    # load fRG data
    DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, Lambda, Delta, U, V, \
    v, Sigma, A, \
    w, K1a, K1p, K1t, \
    w2, v2, K2a, K2p, K2t \
        = load_data(U_NRG, Delta_NRG, filenames)
    
    fs = 18 # font size
    
    lw = 1. #linewidth
    #Plot spectral funciton
    if(typ[0]=="A"):        
        fig, ax = plt.subplots(figsize=(18, 6))

        for i in range(len(filenames)):
            if (filenames[i][-43:-37] == 'static'):
                ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
            else:
                ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:    
            ax.plot(w_NRG/Delta_NRG, Aimp*np.pi*Delta_NRG, 'k:', label='NRG')    
            
        ax.hlines(1., -0.5, 0.5, linestyles='--')
        ax.set_ylabel(r'$(\mathcal{A}(\nu)\pi\Delta$', fontsize=fs)
        
        ax.set_xlim(-5, 5)
        ax.set_xlabel(r'$\nu/\Delta$', fontsize=fs)
        ax.tick_params(axis='both', labelsize=fs-1)
    
        ax.legend(fontsize=fs)
    
    #Plot self energy, component iK
    elif(typ[0:4] == "self"):
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
        
        for i in range(len(filenames)):
            if (filenames[i][-43:-37] == 'static'):                
                ax[0].plot(v[i]/Delta[i], Sigma[i][iK].real/Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                ax[1].plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
            else:
                ax[0].plot(v[i]/Delta[i], Sigma[i][iK].real/Delta[i], '.-', label=labels[i], linewidth=lw)
                ax[1].plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:
            Re_X_NRG = (SE_re - U_NRG / 2) / Delta_NRG
            Im_X_NRG = SE_im / Delta_NRG
            if iK == 1:
                Re_X_NRG = np.zeros(len(w_NRG/Delta_NRG))
                Im_X_NRG = generate_Keldysh_component( Im_X_NRG, w_NRG / Delta_NRG, T[i], stat=-1)
            
            ax[0].plot(w_NRG / Delta_NRG, Re_X_NRG, 'k:', label='NRG')
            ax[1].plot(w_NRG / Delta_NRG, Im_X_NRG, 'k:', label='NRG')
            
            
        if inset:
            axins = zoomed_inset_axes(ax[1], 10, loc=4)
            axins.aspect='equal'
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):                
                    axins.plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    axins.plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw)
            axins.set_xlim(-0.4, 0.4)
            try:
                axins.set_ylim(U_to_ins_xlim[U_NRG], 0.001)
            except:
                axins.set_ylim(-U_NRG/10, 0.001)
            axins.set_xticks(np.linspace(-0.4, 0.4, 5, endpoint=True))
            axins.set_yticks(np.linspace(-0.001, 0.001, 3, endpoint=True))
            plt.xticks(visible=False)  # Not present ticks
            plt.yticks(visible=False)
            
        if(iK==0):
            ax[1].hlines(0, -2, 2, linestyles='--')
            if plot_NRG and inset:
                    axins.plot(w_NRG / Delta_NRG, SE_im / Delta_NRG, 'k:', label='NRG')
            if inset:
                axins.hlines(0, -2, 2, linestyles='--')

            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^R-\Sigma^H)/\Delta$', fontsize=fs)
            ax[1].set_ylabel('$(\mathrm{Im}\Sigma^R)/\Delta$', fontsize=fs)
        else:
            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^K)/\Delta$', fontsize=fs)
            ax[1].set_ylabel('$(\mathrm{Im}\Sigma^K)/\Delta$', fontsize=fs)

            
        for i in range(2):
            ax[i].set_xlim(-20, 20)
            ax[i].set_xlabel(r'$\nu/\Delta$', fontsize=fs)
            ax[i].tick_params(axis='both', labelsize=fs-1)

        ax[0].legend(fontsize=fs)
        #fig.tight_layout(pad=3.0)
    
    #Plot susceptibilities
    elif(typ[0:4] == "susc"):    #susceptibility
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
        sup_indices = ['R', 'K']
        pfs = [-1, +1]
        if typ[6:7]=="sp":      #spin susceptibility
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):
                    ax[0].plot(w[i] / Delta[i],  K1a[i][iK].real / Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*K1a[i][iK].imag / Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    ax[0].plot(w[i] / Delta[i],  K1a[i][iK].real / Delta[i], '.-', label=labels[i], linewidth=lw)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*K1a[i][iK].imag / Delta[i], '.-', label=labels[i], linewidth=lw)
                    
        elif typ[6:7]=="ch":        #charge susceptibility
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):
                    ax[0].plot(w[i] / Delta[i],  (2*K1t[i][0]-K1a[i][iK]).real / Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*(2*K1t[i][0]-K1a[i][iK]).imag / Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    ax[0].plot(w[i] / Delta[i],  (2*K1t[i][0]-K1a[i][iK]).real / Delta[i], '.-', label=labels[i], linewidth=lw) 
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*(2*K1t[i][0]-K1a[i][iK]).imag / Delta[i], '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:
            if typ[6:7] == "sp":
                Re_X_NRG = K1a_re * (U_NRG/Delta_NRG)**2 / 2
                Im_X_NRG = K1a_im * (U_NRG/Delta_NRG)**2 / 2
                
            elif typ[6:7] == "ch":
                Re_X_NRG = (2*K1t_re - K1a_re) * (U_NRG/Delta_NRG)**2 / 2
                Im_X_NRG = (2*K1t_im - K1a_im) * (U_NRG/Delta_NRG)**2 / 2
                
            if iK == 1:
                Re_X_NRG = np.zeros(len(w_NRG/Delta_NRG))
                Im_X_NRG = generate_Keldysh_component( Im_X_NRG, w_NRG / Delta_NRG, T[i], stat=+1)
            
            ax[0].plot(w_NRG / Delta_NRG, Re_X_NRG, 'k:', label='NRG')
            ax[1].plot(w_NRG / Delta_NRG, Im_X_NRG, 'k:', label='NRG')
            
        
        ax[0].set_ylabel(r'$\mathrm{Re}\chi^'+ sup_indices[iK]+'_\mathrm{'+typ[6:7]+'}/\Delta$', fontsize=fs)
        ax[1].set_ylabel(r'$\mathrm{Im}\chi^'+ sup_indices[iK]+'_\mathrm{'+typ[6:7]+'}/\Delta$', fontsize=fs)
        

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs)
            if typ[6:7] == 'sp':
                ax[i].set_xlim(-5,5)
            elif typ[6:7] == 'ch':
                ax[i].set_xlim(10,10)
            ax[i].tick_params(axis='both', labelsize=fs -1)
        ax[0].legend(fontsize=fs)
        fig.tight_layout(pad=3.0)            
    return fig, ax
