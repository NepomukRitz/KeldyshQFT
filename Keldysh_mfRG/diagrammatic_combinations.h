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

#include "vertex.h" // vertex class

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
auto left_same_bare (const Vertex<Q>& vertex, int i1, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2b;
    gamma0 = vertex[1].irred.val(i1, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[1].avertex.K1_valsmooth(i1, w, i_in, vertex[1].tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[1].avertex.K2b_valsmooth(i1, w, vpp, i_in, vertex[1].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[1].pvertex.K1_valsmooth(i1, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[1].pvertex.K2b_valsmooth(i1, w, vpp, i_in);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[1].tvertex.K1_valsmooth(i1, w, i_in, vertex[1].avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[1].tvertex.K2b_valsmooth(i1, w, vpp, i_in, vertex[1].avertex);
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
auto right_same_bare (const Vertex<Q>& vertex, int i3, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2;
    gamma0 = vertex[1].irred.val(i3, i_in);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[1].avertex.K1_valsmooth(i3, w, i_in, vertex[1].tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[1].avertex.K2_valsmooth(i3, w, vpp, i_in, vertex[1].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[1].pvertex.K1_valsmooth(i3, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[1].pvertex.K2_valsmooth(i3, w, vpp, i_in);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[1].tvertex.K1_valsmooth(i3, w, i_in, vertex[1].avertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[1].tvertex.K2_valsmooth(i3, w, vpp, i_in, vertex[1].avertex);
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
auto left_diff_bare (const Vertex<Q>& vertex, int i1, double w, double v, double vpp, int i_in, char channel) -> Q {

    Q K2, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex[1].gammaRb(i1, w, v, vpp, i_in, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex[1].avertex.K2_valsmooth(i1, w, v, i_in, vertex[1].tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[1].avertex.K3_valsmooth(i1, w, v, vpp, i_in, vertex[1].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex[1].pvertex.K2_valsmooth(i1, w, v, i_in);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[1].pvertex.K3_valsmooth(i1, w, v, vpp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex[1].tvertex.K2_valsmooth(i1, w, v, i_in, vertex[1].avertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[1].tvertex.K3_valsmooth(i1, w, v, vpp, i_in, vertex[1].avertex);
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
auto right_diff_bare (const Vertex<Q>& vertex, int i3, double w, double vp, double vpp, int i_in, char channel) -> Q {

    Q K2b, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex[1].gammaRb(i3, w, vpp, vp, i_in, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex[1].avertex.K2b_valsmooth(i3, w, vp, i_in, vertex[1].tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[1].avertex.K3_valsmooth(i3, w, vpp, vp, i_in, vertex[1].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex[1].pvertex.K2b_valsmooth(i3, w, vp, i_in);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[1].pvertex.K3_valsmooth(i3, w, vpp, vp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex[1].tvertex.K2b_valsmooth(i3, w, vp, i_in, vertex[1].avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[1].tvertex.K3_valsmooth(i3, w, vpp, vp, i_in, vertex[1].avertex);
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
auto left_same_bare (const Vertex<Q>& vertex, int i1, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2b;
    gamma0 = vertex[0].irred.val(i1, i_in, spin);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[0].avertex.K1_valsmooth(i1, w, i_in, spin, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].avertex.K2b_valsmooth(i1, w, vpp, i_in, spin, vertex[0].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[0].pvertex.K1_valsmooth(i1, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].pvertex.K2b_valsmooth(i1, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[0].tvertex.K1_valsmooth(i1, w, i_in, spin, vertex[0].avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].tvertex.K2b_valsmooth(i1, w, vpp, i_in, spin, vertex[0].avertex);
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
auto right_same_bare (const Vertex<Q>& vertex, int i3, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2;
    gamma0 = vertex[0].irred.val(i3, i_in, spin);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[0].avertex.K1_valsmooth(i3, w, i_in, spin, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].avertex.K2_valsmooth(i3, w, vpp, i_in, spin, vertex[0].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[0].pvertex.K1_valsmooth(i3, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].pvertex.K2_valsmooth(i3, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[0].tvertex.K1_valsmooth(i3, w, i_in, spin, vertex[0].avertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].tvertex.K2_valsmooth(i3, w, vpp, i_in, spin, vertex[0].avertex);
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
auto left_diff_bare (const Vertex<Q>& vertex, int i1, double w, double v, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex[0].gammaRb(i1, w, v, vpp, i_in, spin, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex[0].avertex.K2_valsmooth(i1, w, v, i_in, spin, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].avertex.K3_valsmooth(i1, w, v, vpp, i_in, spin, vertex[0].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex[0].pvertex.K2_valsmooth(i1, w, v, i_in, spin);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].pvertex.K3_valsmooth(i1, w, v, vpp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex[0].tvertex.K2_valsmooth(i1, w, v, i_in, spin, vertex[0].avertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].tvertex.K3_valsmooth(i1, w, v, vpp, i_in, spin, vertex[0].avertex);
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
auto right_diff_bare (const Vertex<Q>& vertex, int i3, double w, double vp, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2b, K3, gammaRb;
#if DIAG_CLASS >= 2
    gammaRb = vertex[0].gammaRb(i3, w, vpp, vp, i_in, spin, channel);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex[0].avertex.K2b_valsmooth(i3, w, vp, i_in, spin, vertex[0].tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].avertex.K3_valsmooth(i3, w, vpp, vp, i_in, spin, vertex[0].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex[0].pvertex.K2b_valsmooth(i3, w, vp, i_in, spin);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].pvertex.K3_valsmooth(i3, w, vpp, vp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex[0].tvertex.K2b_valsmooth(i3, w, vp, i_in, spin, vertex[0].avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].tvertex.K3_valsmooth(i3, w, vpp, vp, i_in, spin, vertex[0].avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}



#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
