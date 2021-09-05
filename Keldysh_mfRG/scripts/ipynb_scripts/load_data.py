import h5py
import numpy as np

def load_parameters(filename):
    """
    Load system parameters from file <filename>
    """
    with h5py.File(filename, 'r') as f:
        parameters = np.array(f["parameters"])
        Gamma      = parameters[2]
        DIAG_CLASS = parameters[3]
        N_LOOPS    = parameters[4]
        T          = parameters[5]
        mu         = parameters[6]
        U          = parameters[7]
        epsilon    = parameters[8]
        V          = parameters[9]
        Lambdas    = np.array(f["lambdas"])
    return DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, U, V, Lambdas

def load_SE_A(filename, iLambda=0):
    """
    Load self-energy and spectral function from hdf5 file.

    Parameters:
        filename: name of the file from which data should be read
        iLambda:  Lambda iteration that should be read off (default: 0)

    Returns:
        v:          vector of frequencies at which result is saved
        SigmaR:     retarded self-energy
        SigmaK:     Keldysh self-energy
        A:          spectral function
    """
    with h5py.File(filename, 'r') as f:
        v = np.array(f["ffreqs"])
        self = np.array(f["selflist"][iLambda])

    _, _, _, epsilon, _, Gamma, U, _, Lambdas = load_parameters(filename)

    nv = len(v)
    Sigma = np.zeros((2, nv), dtype=complex)
    for iK in range(2):
        for iv in range(len(v)):
            Sigma[iK][iv] = self[iK * nv + iv][0] + 1j * self[iK * nv + iv][1]

    G = 1/(v - epsilon + (Gamma+Lambdas[iLambda])/2*1j - Sigma[0])
    A = -np.imag(G)/np.pi

    for iv in range(len(v)):
        Sigma[0][iv] -= U/2.
    
    return v, Sigma, A

def load_K1(filename, iLambda=0):
    """
        Load K1 class from hdf5 file.

        Parameters:
            filename: name of the file from which data should be read
            iLambda:  Lambda iteration that should be read off (default: 0)

        Returns:
            w:   vector of frequencies at which result is saved
            K1a: K1 in the a channel
            K1p: K1 in the p channel
            K1t: K1 in the t channel
        """
    with h5py.File(filename, 'r') as f:
        w = np.array(f["bfreqs"])
        K1_a = np.array(f["K1_a"][iLambda])
        K1_p = np.array(f["K1_p"][iLambda])
        K1_t = np.array(f["K1_t"][iLambda])

    nw = len(w)
    K1a = np.zeros((2, nw), dtype=complex)
    K1p = np.zeros((2, nw), dtype=complex)
    K1t = np.zeros((2, nw), dtype=complex)
    for iK in range(2):
        for iw in range(nw):
            K1a[iK][iw] = K1_a[iK * nw + iw][0] + 1j * K1_a[iK * nw + iw][1]
            K1p[iK][iw] = K1_p[iK * nw + iw][0] + 1j * K1_p[iK * nw + iw][1]
            K1t[iK][iw] = K1_t[iK * nw + iw][0] + 1j * K1_t[iK * nw + iw][1]

    return w, K1a, K1p, K1t

def load_K2(filename, iLambda=0):
    """
    Same as load_K1
    """
    with h5py.File(filename, 'r') as f:
        w2 = np.array(f["bfreqs2"])
        v2 = np.array(f["ffreqs2"])
        K2_a = np.array(f["K2_a"][iLambda])
        K2_p = np.array(f["K2_p"][iLambda])
        K2_t = np.array(f["K2_t"][iLambda])
    nw2 = len(w2)
    nv2 = len(v2)
    K2a = np.zeros((5, nw2, nv2), dtype=complex)
    K2p = np.zeros((5, nw2, nv2), dtype=complex)
    K2t = np.zeros((5, nw2, nv2), dtype=complex)

    for iK in range(5):
        for iw in range(nv2):
            for iv in range(nv2):
                K2a[iK][iw][iv] = K2_a[iK * nw2 * nv2 + iw * nv2 + iv][0] + 1j * K2_a[iK * nw2 * nv2 + iw * nv2 + iv][1]
                K2p[iK][iw][iv] = K2_p[iK * nw2 * nv2 + iw * nv2 + iv][0] + 1j * K2_p[iK * nw2 * nv2 + iw * nv2 + iv][1]
                K2t[iK][iw][iv] = K2_t[iK * nw2 * nv2 + iw * nv2 + iv][0] + 1j * K2_t[iK * nw2 * nv2 + iw * nv2 + iv][1]
    return w2, v2, K2a, K2p, K2t

def load_data(U_NRG, Delta_NRG, filenames):
    """
    Load all fRG data
    
    Parameters:
        U_NRG:     value of the interaction U for NRG data at which to load fRG data
        Delta_NRG: effective hybridization of NRG data (usually Gamma_NRG=1 ==> Delta_NRG=0.5)
        filenames: list of file names from which the data should be obtained
    """
    
    DIAG_CLASS_list, N_LOOPS_list = [], []
    T_list, epsilon_list, mu_list = [], [], []
    Gamma_list, Lambda_list, Delta_list, U_list, V_list = [], [], [], [], []
    v_list, Sigma_list, A_list = [], [], []
    w_list, K1a_list, K1p_list, K1t_list = [], [], [], []
    w2_list, v2_list, K2a_list, K2p_list, K2t_list = [], [], [], [], []

    for filename in filenames:
        # load parameters
        DIAG_CLASS, N_LOOPS, T, epsilon, mu, Gamma, U, V, Lambdas = load_parameters(filename)

        # determine effective hybridization Delta
        Deltas = (Gamma + Lambdas) / 2

        # check at which Lambda iteration U/Delta matches between fRG and the desired NRG value
        if(U_NRG == 0):
            iLambda = 0
        else:
            try:
                iLambda = np.where(U/Deltas == U_NRG/Delta_NRG)[0][0]
            except IndexError:
                iLambda = np.argmin(np.abs(U/Deltas - U_NRG/Delta_NRG))
                print("Actual U/Delta-value = " + str( (2*U)/(Lambdas[iLambda] + Gamma)))

        # fix Lambda and effective hybridization Delta
        Lambda = Lambdas[iLambda]
        Delta = (Gamma + Lambda) / 2

        # load self-energy and spectral function
        v, Sigma, A = load_SE_A(filename, iLambda=iLambda)

        # load K1
        w, K1a, K1p, K1t = load_K1(filename, iLambda=iLambda)

        # load K2
        if DIAG_CLASS >= 2:
            w2, v2, K2a, K2p, K2t = load_K2(filename, iLambda=iLambda)
        else:
            w2, v2, K2a, K2p, K2t = [], [], [], [], []

        # add all data to lists
        DIAG_CLASS_list.append(DIAG_CLASS)
        N_LOOPS_list.append(N_LOOPS)
        T_list.append(T)
        epsilon_list.append(epsilon)
        mu_list.append(mu)
        Gamma_list.append(Gamma)
        Lambda_list.append(Lambda)
        Delta_list.append(Delta)
        U_list.append(U)
        V_list.append(V)

        v_list.append(v)
        Sigma_list.append(Sigma)
        A_list.append(A)
        w_list.append(w)
        K1a_list.append(K1a)
        K1p_list.append(K1p)
        K1t_list.append(K1t)
        w2_list.append(w2)
        v2_list.append(v2)
        K2a_list.append(K2a)
        K2p_list.append(K2p)
        K2t_list.append(K2t)

    return DIAG_CLASS_list, N_LOOPS_list, \
        T_list, epsilon_list, mu_list, Gamma_list, Lambda_list, Delta_list, U_list, V_list, \
        v_list, Sigma_list, A_list, \
        w_list, K1a_list, K1p_list, K1t_list, \
        w2_list, v2_list, K2a_list, K2p_list, K2t_list


def load_SIAM_NRG(U, path_NRG):
    """
    Load NRG reference data
    """
    with h5py.File(path_NRG + "SIAM_NRG.h5", 'r') as f:
        Gamma_NRG = np.array(f['Gamma'])
        w_NRG  = np.array(f["w"])
        Aimp   = np.array(f["U=" + str(U)]["Aimp"])
        SE_re  = np.array(f["U=" + str(U)]["SE_re"])
        SE_im  = np.array(f["U=" + str(U)]["SE_im"])
        K1a_re = np.array(f["U=" + str(U)]["K1a_div_by_U_squared_re"])
        K1a_im = np.array(f["U=" + str(U)]["K1a_div_by_U_squared_im"])
        K1p_re = np.array(f["U=" + str(U)]["K1p_div_by_U_squared_re"])
        K1p_im = np.array(f["U=" + str(U)]["K1p_div_by_U_squared_im"])
        K1t_re = np.array(f["U=" + str(U)]["K1t_div_by_U_squared_re"])
        K1t_im = np.array(f["U=" + str(U)]["K1t_div_by_U_squared_im"])

    return Gamma_NRG, w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im

def load_SIAM_NRG_finite_T(U, T, path_NRG):
    """
    Load NRG reference data
    """
    if T==0.01:
        return load_SIAM_NRG(U, path_NRG)
    
    with h5py.File(path_NRG + "SIAM_NRG_finite_T.h5", 'r') as f:
        Gamma_NRG = np.array(f['Gamma'])
        w_NRG  = np.array(f["w"])
        Aimp   = np.array(f["U=" + str(U)+'_T=' + str(T)]["Aimp"])
        SE_re  = np.array(f["U=" + str(U)+'_T=' + str(T)]["SE_re"])
        SE_im  = np.array(f["U=" + str(U)+'_T=' + str(T)]["SE_im"])
        K1a_re = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1a_div_by_U_squared_re"])
        K1a_im = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1a_div_by_U_squared_im"])
        K1p_re = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1p_div_by_U_squared_re"])
        K1p_im = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1p_div_by_U_squared_im"])
        K1t_re = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1t_div_by_U_squared_re"])
        K1t_im = np.array(f["U=" + str(U)+'_T=' + str(T)]["K1t_div_by_U_squared_im"])

    return Gamma_NRG, w_NRG, Aimp, SE_re, SE_im, K1a_re, K1a_im, K1p_re, K1p_im, K1t_re, K1t_im