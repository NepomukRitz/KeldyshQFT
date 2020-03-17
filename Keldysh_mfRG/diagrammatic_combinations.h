/**
 * Collect combinations of diagrammatic classes (Gamma0, K1, K2, K3) that contribute to the left or right vertex
 * in bubbles in the flow equations for K1, K2, K3.
 * There are only two combinations that occur for both the left and the right vertex:
 *  - Those diag. classes that connect to the same bare vertex on the left/right side, respectively. These are
 *      Gamma0, K1, K2b on the left (left_same_bare) and
 *      Gamma0, K1, K2 on the right (right_same_bare).
 *  - Those diag. classes that connect to different bare vertices on the left/right side, respectively. These are
 *      K2, K3, gamma_bar{r} on the left (left_diff_bare) and
 *      K2b, K_3, gamma_bar{r} on the right (right_diff_bare).
 */

#ifndef KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
#define KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H

#include "vertex.h"
// TODO: vval -> val, vvalsmooth -> valsmooth

/**
 * Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
 * (for spin component densvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto left_same_bare (const Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2b;
    gamma0 = vertex.densvertex.irred.vval(i1, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.avertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.pvertex.K1_vvalsmooth(i1, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.tvertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return gamma0 + K1 + K2b;
}

/**
 * Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
 * (for spin component densvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i3      : Keldysh index for this (right) vertex
 * @param w       : bosonic frequency omega
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto right_same_bare (const Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2; // TODO: check frequency
    gamma0 = vertex.densvertex.irred.vval(i3, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.avertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.avertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.pvertex.K1_vvalsmooth(i3, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.pvertex.K2_vvalsmooth(i3, w, vpp, i_in);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.tvertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.tvertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return gamma0 + K1 + K2;
}

/**
 * Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
 * (for spin component densvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param v       : fermionic frequency nu
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto left_diff_bare (const Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in, char channel) -> Q {

    Q K2, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gammaRb;
}

/**
 * Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
 * (for spin component densvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param vp      : fermionic frequency nu'
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto right_diff_bare (const Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in, char channel) -> Q {

    Q K2b, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}


/**
 * Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
 * (for spin component spinvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param spin    : spin component: 0 = V, 1 = hat{V}
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return        : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto left_same_bare (const Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2b;
    gamma0 = vertex.spinvertex.irred.vval(i1, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.avertex.K1_vvalsmooth(i1, w, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.pvertex.K1_vvalsmooth(i1, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.tvertex.K1_vvalsmooth(i1, w, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return gamma0 + K1 + K2b;
}

/**
 * Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
 * (for spin component spinvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i3      : Keldysh index for this (right) vertex
 * @param w       : bosonic frequency omega
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param spin    : spin component: 0 = V, 1 = hat{V}
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return        : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto right_same_bare (const Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2;
    gamma0 = vertex.spinvertex.irred.vval(i3, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.avertex.K1_vvalsmooth(i3, w, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.avertex.K2_vvalsmooth(i3, w, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.pvertex.K1_vvalsmooth(i3, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.pvertex.K2_vvalsmooth(i3, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.tvertex.K1_vvalsmooth(i3, w, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.tvertex.K2_vvalsmooth(i3, w, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return gamma0 + K1 + K2;
}

/**
 * Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
 * (for spin component spinvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param v       : fermionic frequency nu
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param spin    : spin component: 0 = V, 1 = hat{V}
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto left_diff_bare (const Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex.spinvertex.gammaRb(i1, w, v, vpp, i_in, spin, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in, spin);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gammaRb;
}

/**
 * Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
 * (for spin component densvertex)
 * @tparam Q      : return type (comp or double)
 * @param vertex  : vertex from which diagrammatic classes are taken
 * @param i1      : Keldysh index for this (left) vertex
 * @param w       : bosonic frequency omega
 * @param vp      : fermionic frequency nu'
 * @param vpp     : fermionic frequency nu'' (to be integrated over)
 * @param i_in    : internal structure index
 * @param channel : diagrammatic channel ('a', 'p', 't')
 * @return Q      : result of the corresponding combination of diag. classes evaluated at the above arguments
 */
template <typename Q>
auto right_diff_bare (const Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2b, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex.spinvertex.gammaRb(i3, w, vpp, vp, i_in, spin, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in, spin);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}



#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
