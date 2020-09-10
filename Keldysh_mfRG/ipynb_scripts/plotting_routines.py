import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from load_data import load_data, load_SIAM_NRG

# use latex and serif fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot_self_energy(U_NRG, Delta_NRG, filenames, labels, path="../cmake-build-debug-remote/"):
    """
    Two-panel plot of real and imaginary part of the retarded self-energy

    Parameters:
        U_NRG:     value of the interaction U at which to plot result
        Delta_NRG: effective hybridization of NRG data (usually Gamma_NRG=1 ==> Delta_NRG=0.5)
        filenames: list of file names from which the data should be plotted
        labels:    list of labels for the legend, specifying the data sets in filenames
        path:      path where to find the files (default: Keldysh_mfRG/cmake-build-debug-remote/)

    Returns:
        figure and axis handle
    """

    # load fRG data
    DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, Lambda, Delta, U, V, \
    v, SigmaR, SigmaK, A, \
    w, K1a, K1p, K1t, \
    w2, v2, K2a, K2p, K2t \
        = load_data(U_NRG, Delta_NRG, filenames)

    # load NRG reference data
    w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im = load_SIAM_NRG(U_NRG)

    fs = 18 # font size
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

    for i in range(len(v)):
        ax[0].plot(v[i]/Delta[i], (SigmaR[i].real-U[i]/2)/Delta[i], '.-', label=labels[i])
        ax[1].plot(v[i]/Delta[i], SigmaR[i].imag/Delta[i], '.-', label=labels[i])

    ax[0].plot(w_NRG / Delta_NRG, (SE_re - U_NRG / 2) / Delta_NRG, 'k:', label='NRG')
    ax[1].plot(w_NRG / Delta_NRG, SE_im / Delta_NRG, 'k:', label='NRG')

    ax[0].set_ylabel('$(\mathrm{Re}\Sigma^R-U/2)/\Delta$', fontsize=fs)
    ax[1].set_ylabel('$(\mathrm{Im}\Sigma^R-U/2)/\Delta$', fontsize=fs)

    for i in range(2):
        ax[i].set_xlabel(r'$\nu/\Delta$', fontsize=fs)
        ax[i].tick_params(axis='both', labelsize=fs-1)

    ax[0].legend(fontsize=fs)
    fig.tight_layout(pad=3.0)

    return fig, ax

def plot_spin_susceptibility(U_NRG, Delta_NRG, filenames, labels, path="../cmake-build-debug-remote/"):
    """
    Two-panel plot of real and imaginary part of the spin susceptibility

    Similar to plot_self_energy
    """

    # load fRG data
    DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, Lambda, Delta, U, V, \
    v, SigmaR, SigmaK, A, \
    w, K1a, K1p, K1t, \
    w2, v2, K2a, K2p, K2t \
        = load_data(U_NRG, Delta_NRG, filenames)

    # load NRG reference data
    w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im = load_SIAM_NRG(U_NRG)

    fs = 18  # font size
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

    for i in range(len(v)):
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
        ax[i].tick_params(axis='both', labelsize=fs - 1)

    ax[0].legend(fontsize=fs)
    fig.tight_layout(pad=3.0)

    return fig, ax
