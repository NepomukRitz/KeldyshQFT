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
        indices = np.where(np.abs(res) > 10)[0]
        for i in indices:
            res[i] = np.sign(res[i])*10
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
        label = r'$\mathcal{K}_' + diag_class + '$'
        if (filename[len(path_fRG)+14 : len(path_fRG)+20] == 'static'):
            label +=' -- sf'
        
        elif(diag_class == str(2)):
            label += ' ' + filename[len(path_fRG)+3] + '-loop'
            
        labels.append(label)
    return labels


def plot_fRG(typ, iK, filenames, labels, NRG_info, inset):
    U_NRG = NRG_info[0]
    path_NRG = NRG_info[1]
    possible_NRG = NRG_info[2]
    
    U_to_ins_ylim  = {0.1: 0.0001, 0.5:0.003, 1:0.008, 1.5: 0.02, 2:0.05, 2.5:0.07, 2.9:0.1, 3:0.1}
    U_str = "%.1f" % (2*U_NRG)
    
    
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
    
    K1_dic = {'K1a': K1a, 'K1p': K1p, 'K1t': K1t}
    
    fs = 26 # font size
    fs_leg = fs - 6
    fs_ticks = fs-2
    fs_ins = fs_leg
    plt.rcParams['legend.title_fontsize'] = fs-4
    
    lw = 3 #linewidth
    #Plot spectral funciton
    if(typ[0]=="A"):        
        fig, ax = plt.subplots(figsize=(15, 5))

        for i in range(len(filenames)):
            if (filenames[i][-43:-37] == 'static'):
                ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
            else:
                ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:    
            ax.plot(w_NRG/Delta_NRG, Aimp*np.pi*Delta_NRG, 'k:', label='NRG')    
            
        ax.hlines(1., -0.5, 0.5, linestyles='--')
        ax.set_ylabel(r'$\mathcal{A}(\nu)\pi\Delta$', fontsize=fs)
        
        ax.set_ylim(-0.1, 1.15)
        ax.set_xlim(-5, 5)
        ax.set_xlabel(r'$\nu/\Delta$', fontsize=fs)
        ax.tick_params(axis='both', labelsize=fs_ticks)
        
                    
        if inset:
            axins = zoomed_inset_axes(ax, 2.1, loc=2, borderpad = 4)
            axins.aspect='equal'
            axins.hlines(1., -0.5, 0.5, linestyles='--')
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):                
                    axins.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    axins.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], '.-', label=labels[i], linewidth=lw)
            if plot_NRG:    
                axins.plot(w_NRG/Delta_NRG, Aimp*np.pi*Delta_NRG, 'k:', label='NRG')                      
                
                
            axins.set_ylim(0.9, 1.1)
            axins.set_yticks(np.linspace(0.9, 1.1, 3, endpoint=True))    
            axins.set_xlim(-0.4, 0.4)    
            axins.set_xticks(np.linspace(-0.4, 0.4, 3, endpoint=True))
            axins.tick_params(axis='both', labelsize=fs_ins)
            plt.xticks(visible=True)  # Present ticks
            plt.yticks(visible=True)
            mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
        ax.legend(fontsize=fs_leg, title=r'$U/\Delta=$' + U_str)
    
    #Plot self energy, component iK
    elif(typ[0:4] == "self"):
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        
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
            axins = zoomed_inset_axes(ax[1], 7, loc=4, borderpad=2.5)
            axins.aspect='equal'
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):                
                    axins.plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    axins.plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], '.-', label=labels[i], linewidth=lw)
            try:
                ylim = U_to_ins_ylim[U_NRG]
            except KeyError:
                ylim = U_NRG/20
                
                
            axins.set_ylim(-ylim, 0.2*ylim)
            axins.set_yticks(np.linspace(-ylim, 0, 2, endpoint=True))    
            axins.set_xlim(-0.5, 0.5)    
            axins.set_xticks(np.linspace(-0.5, 0.5, 3, endpoint=True))
            axins.tick_params(axis='both', labelsize=fs_ins-4)

            
            plt.xticks(visible=True)  # Present ticks
            plt.yticks(visible=True)
            
        if(iK==0):
            ax[1].hlines(0, -2, 2, linestyles='--')
            if plot_NRG and inset:
                    axins.plot(w_NRG / Delta_NRG, SE_im / Delta_NRG, 'k:', label='NRG')
            if inset:
                axins.hlines(0, -2, 2, linestyles='--')

            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^R-\Sigma^H)/\Delta$', fontsize=fs+2)
            ax[1].set_ylabel('$\mathrm{Im}\Sigma^R/\Delta$', fontsize=fs+2)
        else:
            ax[0].set_ylabel('$\mathrm{Re}\Sigma^K/\Delta$', fontsize=fs+2)
            ax[1].set_ylabel('$\mathrm{Im}\Sigma^K/\Delta$', fontsize=fs+2)

            
        for i in range(2):
            ax[i].set_xlim(-20, 20)
            ax[i].set_xlabel(r'$\nu/\Delta$', fontsize=fs+2)
            ax[i].tick_params(axis='both', labelsize=fs_ticks)

        ax[0].legend(fontsize=fs_leg, loc='upper left', title=r'$U/\Delta=$' + U_str)
        plt.tight_layout()
        
    #Plot susceptibilities
    elif(typ[0:4] == "susc"):    #susceptibility
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        sup_indices = ['R', 'K']
        pfs = [-1, +1]
        susc = typ[5:7]
        if susc=="sp":      #spin susceptibility
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):
                    ax[0].plot(w[i] / Delta[i],  K1a[i][iK].real *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*K1a[i][iK].imag *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    ax[0].plot(w[i] / Delta[i],  K1a[i][iK].real *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*K1a[i][iK].imag *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw)
                    
        elif susc=="ch":        #charge susceptibility
            for i in range(len(filenames)):
                if (filenames[i][-43:-37] == 'static'):
                    ax[0].plot(w[i] / Delta[i], (2*K1t[i][iK]-K1a[i][iK]).real *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*(2*K1t[i][iK]-K1a[i][iK]).imag *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw*0.5)
                else:
                    ax[0].plot(w[i] / Delta[i], (2*K1t[i][iK]-K1a[i][iK]).real *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw) 
                    ax[1].plot(w[i] / Delta[i], pfs[iK]*(2*K1t[i][iK]-K1a[i][iK]).imag *np.pi* Delta[i], '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:
            if susc == "sp":
                Re_X_NRG = K1a_re 
                Im_X_NRG = K1a_im 
                
            elif susc == "ch":
                Re_X_NRG = (2*K1t_re - K1a_re)
                Im_X_NRG = (2*K1t_im - K1a_im)
                
            if iK == 1:
                Re_X_NRG = np.zeros(len(w_NRG/Delta_NRG))
                Im_X_NRG = generate_Keldysh_component( Im_X_NRG, w_NRG, T[i], stat=+1)
            
            ax[0].plot(w_NRG / Delta_NRG, Re_X_NRG * 2*U_NRG**2* np.pi*Delta[i]**2, 'k:', label='NRG') 
            ax[1].plot(w_NRG / Delta_NRG, Im_X_NRG * 2*U_NRG**2* np.pi*Delta[i]**2, 'k:', label='NRG')
            
        
        ax[0].set_ylabel(r'$\mathrm{Re}\chi^'+ sup_indices[iK]+'_\mathrm{'+susc+'}\pi\Delta$', fontsize=fs+2)
        ax[1].set_ylabel(r'$\mathrm{Im}\chi^'+ sup_indices[iK]+'_\mathrm{'+susc+'}\pi\Delta$', fontsize=fs+2)
        

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs+2)
            if susc == 'sp':
                ax[i].set_xlim(-5,5)
            elif susc == 'ch':
                ax[i].set_xlim(-10,10)
            ax[i].tick_params(axis='both', labelsize=fs_ticks)
        ax[0].legend(fontsize=fs_leg, title=r'$U/\Delta=$' + U_str)
        fig.tight_layout()  
        
        
    #Plot vertex K1
    elif(typ[0:2] == "K1"):    
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        sup_indices = ['R', 'K']
        r = typ[2]
        pfs = [-1, +1]
        K1 = K1_dic[typ]
        for i in range(len(filenames)):
            if (filenames[i][-43:-37] == 'static'):
                ax[0].plot(w[i] / Delta[i], K1[i][iK].real, '.-', label=labels[i], linewidth=lw*0.5)
                ax[1].plot(w[i] / Delta[i], pfs[iK]*K1[i][iK].imag, '.-', label=labels[i], linewidth=lw*0.5)
            else:
                ax[0].plot(w[i] / Delta[i], K1[i][iK].real, '.-', label=labels[i], linewidth=lw)
                ax[1].plot(w[i] / Delta[i], pfs[iK]*K1[i][iK].imag, '.-', label=labels[i], linewidth=lw)
            
        if plot_NRG:
            if r == "a":
                Re_X_NRG = K1a_re 
                Im_X_NRG = K1a_im 
                
            elif r == "p":
                Re_X_NRG = K1p_re 
                Im_X_NRG = K1p_im 
                
            elif r == "t":
                Re_X_NRG = K1t_re 
                Im_X_NRG = K1t_im 
                
            if iK == 1:
                Re_X_NRG = np.zeros(len(w_NRG/Delta_NRG))
                Im_X_NRG = generate_Keldysh_component(Im_X_NRG, w_NRG, T[i], stat=+1)
            
            ax[0].plot(w_NRG / Delta_NRG, Re_X_NRG * (U_NRG)**2 / 2, 'k:', label='NRG')
            ax[1].plot(w_NRG / Delta_NRG, Im_X_NRG * (U_NRG)**2 / 2, 'k:', label='NRG')
            
        
        ax[0].set_ylabel(r'$\mathrm{Re}{\mathcal{K}_1^'+ sup_indices[iK]+'}_\mathrm{'+r+'}$', fontsize=fs)
        ax[1].set_ylabel(r'$\mathrm{Im}{\mathcal{K}_1^'+ sup_indices[iK]+'}_\mathrm{'+r+'}$', fontsize=fs)
        

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs)
            if r == 'a' or 't':
                ax[i].set_xlim(-5,5)
            elif r == 'p':
                ax[i].set_xlim(-10,10)
            ax[i].tick_params(axis='both', labelsize=fs_ticks)
        ax[0].legend(fontsize=fs_leg, title=r'$U/\Delta=$' + U_str)
        fig.tight_layout() 
        
    return fig, ax




def plot_fRG_V_runs(typ, iK, filenames, U_NRG, labels, run, kind):
    if (typ == "self" and iK ==0) or (typ=='A'):
        inset = True
    else:
        inset = False
    
    U_to_ins_ylim  = {0.1: 0.0001, 0.5:0.003, 1:0.008, 1.5: 0.02, 2:0.05, 2.5:0.07, 2.9:0.1, 3:0.1}
    U_str = "%.1f" % (2*U_NRG)

    
    color_dic = {'K1' : 'tab:blue', 'K1_sf' : 'tab:orange', 'K2_1L' : 'tab:green', 'K2_2L' : 'tab:red'}
    title_dic = {'K1' : r'$\mathcal{K}_1$ -- ' + kind,\
                 'K1_sf' : r'$\mathcal{K}_1$ -sf -- ' + kind,\
                 'K2_1L' : r'$\mathcal{K}_2$ 1-L -- ' + kind,\
                 'K2_2L' : r'$\mathcal{K}_2$ 2-loop -- '+'$U/\Delta= '+U_str+'$\n $\quad$'+ kind}
    c = color_dic[run]
    
    styles = ['solid', 'dashed', 'dashdot', 'dotted', ' ', '']
    
    Delta_NRG = 0.5
    
    # load fRG data
    DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, Lambda, Delta, U, V, \
    v, Sigma, A, \
    w, K1a, K1p, K1t, \
    w2, v2, K2a, K2p, K2t \
        = load_data(U_NRG, Delta_NRG, filenames)
    
    K1_dic = {'K1a': K1a, 'K1p': K1p, 'K1t': K1t}
    
    fs = 26 # font size
    fs_leg = fs - 2
    fs_ticks = fs-2
    fs_ins = fs_leg
    
    lw=3
    
    plt.rcParams['legend.title_fontsize'] = fs-2
    
    #Plot spectral funciton
    if(typ[0]=="A"):        
        fig, ax = plt.subplots(figsize=(15, 5))

        for i in range(len(filenames)):
            ax.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], linestyle=styles[i], label=labels[i], color=c, linewidth=lw)
            
        ax.hlines(1., -0.5, 0.5, linestyles='--')
        ax.set_ylabel(r'$\mathcal{A}(\nu)\pi\Delta$', fontsize=fs)
        
        ax.set_ylim(-0.1, 1.15)
        ax.set_xlim(-5, 5)
        ax.set_xlabel(r'$\nu/\Delta$', fontsize=fs)
        ax.tick_params(axis='both', labelsize=fs_ticks)
        
                    
        if inset:
            axins = zoomed_inset_axes(ax, 2.1, loc=2, borderpad = 4)
            axins.aspect='equal'
            axins.hlines(1., -0.5, 0.5, linestyles='--')
            for i in range(len(filenames)):
                axins.plot(v[i]/Delta[i], A[i]*np.pi*Delta[i], linestyle=styles[i], label=labels[i], color=c, linewidth=lw)
                
            axins.set_ylim(0.9, 1.1)
            axins.set_yticks(np.linspace(0.9, 1.1, 3, endpoint=True))    
            axins.set_xlim(-0.4, 0.4)    
            axins.set_xticks(np.linspace(-0.4, 0.4, 3, endpoint=True))
            axins.tick_params(axis='both', labelsize=fs_ins)
            plt.xticks(visible=True)  # Present ticks
            plt.yticks(visible=True)
            mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
        l = ax.legend(fontsize=fs_leg, title=title_dic[run])
        plt.setp(l.get_title(), multialignment='center')
    
    #Plot self energy, component iK
    elif(typ[0:4] == "self"):
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
        
        for i in range(len(filenames)):
            ax[0].plot(v[i]/Delta[i], Sigma[i][iK].real/Delta[i], linestyle=styles[i], label=labels[i], color=c, linewidth=lw)
            ax[1].plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], linestyle=styles[i], label=labels[i], color=c, linewidth=lw)            
            
        if inset:
            axins = zoomed_inset_axes(ax[1], 7, loc=4, borderpad=2.9)
            axins.aspect='equal'
            for i in range(len(filenames)):
                axins.plot(v[i]/Delta[i], Sigma[i][iK].imag/Delta[i], linestyle=styles[i], label=labels[i], color=c, linewidth=lw)
            try:
                ylim = U_to_ins_ylim[U_NRG]
            except KeyError:
                ylim = U_NRG/20
                
                
            axins.set_ylim(-ylim, 0.2*ylim)
            axins.set_yticks(np.linspace(-ylim, 0, 2, endpoint=True))    
            axins.set_xlim(-0.5, 0.5)    
            axins.set_xticks(np.linspace(-0.5, 0.5, 3, endpoint=True))
            axins.tick_params(axis='both', labelsize=fs_ticks-4)

            
            plt.xticks(visible=True)  # Present ticks
            plt.yticks(visible=True)
            
        if(iK==0):
            ax[1].hlines(0, -2, 2, linestyles='--')
            if inset:
                axins.hlines(0, -2, 2, linestyles='--')

            ax[0].set_ylabel('$(\mathrm{Re}\Sigma^R-\Sigma^H)/\Delta$', fontsize=fs+2)
            ax[1].set_ylabel('$\mathrm{Im}\Sigma^R/\Delta$', fontsize=fs+2)
        else:
            ax[0].set_ylabel('$\mathrm{Re}\Sigma^K/\Delta$', fontsize=fs+2)
            ax[1].set_ylabel('$\mathrm{Im}\Sigma^K/\Delta$', fontsize=fs+2)

            
        for i in range(2):
            ax[i].set_xlim(-20, 20)
            ax[i].set_xlabel(r'$\nu/\Delta$', fontsize=fs+2)
            ax[i].tick_params(axis='both', labelsize=fs_ticks)

        l = ax[0].legend(title=title_dic[run], fontsize=fs_leg, loc='upper left')
        plt.setp(l.get_title(), multialignment='center')
        fig.tight_layout()
    
    #Plot susceptibilities
    elif(typ[0:4] == "susc"):    #susceptibility
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
        sup_indices = ['R', 'K']
        pfs = [-1, +1]
        susc = typ[5:7]
        if susc=="sp":      #spin susceptibility
            for i in range(len(filenames)):
                ax[0].plot(w[i] / Delta[i],  K1a[i][iK].real * np.pi* Delta[i],  linestyle=styles[i], label=labels[i], color=c)
                ax[1].plot(w[i] / Delta[i], pfs[iK]*K1a[i][iK].imag * np.pi* Delta[i],  linestyle=styles[i], label=labels[i], color=c)
                    
        elif susc=="ch":        #charge susceptibility
            for i in range(len(filenames)):
                ax[0].plot(w[i] / Delta[i], (2*K1t[i][iK]-K1a[i][iK]).real * np.pi* Delta[i],  linestyle=styles[i], label=labels[i], color=c)
                ax[1].plot(w[i] / Delta[i], pfs[iK]*(2*K1t[i][iK]-K1a[i][iK]).imag * np.pi* Delta[i],  linestyle=styles[i], label=labels[i], color=c)
            
        ax[0].set_ylabel(r'$\mathrm{Re}\chi^'+ sup_indices[iK]+'_\mathrm{'+susc+'}\pi\Delta$', fontsize=fs)
        ax[1].set_ylabel(r'$\mathrm{Im}\chi^'+ sup_indices[iK]+'_\mathrm{'+susc+'}\pi\Delta$', fontsize=fs)
        

        for i in range(2):
            ax[i].set_xlabel(r'$\omega/\Delta$', fontsize=fs)
            if susc == 'sp':
                ax[i].set_xlim(-5,5)
            elif susc == 'ch':
                ax[i].set_xlim(-10,10)
            ax[i].tick_params(axis='both', labelsize=fs -1)
        l = ax[0].legend(fontsize=fs, title=title_dic[run])
        plt.setp(l.get_title(), multialignment='center')
        fig.tight_layout(pad=3.0)  
        
    return fig, ax
