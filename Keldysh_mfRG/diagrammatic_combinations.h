// ATTENTION: DEPRECATED! FULLY MOVED INTO vertex.h
// TODO: Remove dependency on this file from correctionFunctions.h, then remove this file.

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
    VertexInput input (i1, w, 0., vpp, i_in, spin, channel);
    Q gamma0, K1, K2b;
    gamma0 = vertex[0].irred.val(i1, i_in, spin);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[0].avertex.valsmooth(k1, input, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].avertex.valsmooth(k2b, input, vertex[0].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[0].pvertex.valsmooth(k1, input);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].pvertex.valsmooth(k2b, input);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[0].tvertex.valsmooth(k1, input, vertex[0].avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex[0].tvertex.valsmooth(k2b, input, vertex[0].avertex);
#endif
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
#if DIAG_CLASS <= 1
    VertexInput input_p = input;
    VertexInput input_at = input;
    input_p.w = 2*glb_mu;
    input_at.w = 0.;

    switch (channel) {
        case 'a':
            K1 += vertex[0].pvertex.valsmooth(k1, input_p)
                  + vertex[0].tvertex.valsmooth(k1, input_at, vertex[0].avertex);
            break;
        case 'p':
            K1 += vertex[0].avertex.valsmooth(k1, input_at, vertex[0].tvertex)
                  + vertex[0].tvertex.valsmooth(k1, input_at, vertex[0].avertex);
            break;
        case 't':
            K1 += vertex[0].avertex.valsmooth(k1, input_at, vertex[0].tvertex)
                  + vertex[0].pvertex.valsmooth(k1, input_p);
            break;
        default: ;
    }
#endif
#endif
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
    VertexInput input (i3, w, vpp, 0., i_in, spin, channel);
    Q gamma0, K1, K2;
    gamma0 = vertex[0].irred.val(i3, i_in, spin);

    switch (channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = vertex[0].avertex.valsmooth(k1, input, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].avertex.valsmooth(k2, input, vertex[0].tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = vertex[0].pvertex.valsmooth(k1, input);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].pvertex.valsmooth(k2, input);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = vertex[0].tvertex.valsmooth(k1, input, vertex[0].avertex);
#endif
#if DIAG_CLASS >=2
            K2 = vertex[0].tvertex.valsmooth(k2, input, vertex[0].avertex);
#endif
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
#if DIAG_CLASS <= 1
    VertexInput input_p = input;
    VertexInput input_at = input;
    input_p.w = 2*glb_mu;
    input_at.w = 0.;

    switch (channel) {
        case 'a':
            K1 += vertex[0].pvertex.valsmooth(k1, input_p)
                  + vertex[0].tvertex.valsmooth(k1, input_at, vertex[0].avertex);
            break;
        case 'p':
            K1 += vertex[0].avertex.valsmooth(k1, input_at, vertex[0].tvertex)
                  + vertex[0].tvertex.valsmooth(k1, input_at, vertex[0].avertex);
            break;
        case 't':
            K1 += vertex[0].avertex.valsmooth(k1, input_at, vertex[0].tvertex)
                  + vertex[0].pvertex.valsmooth(k1, input_p);
            break;
        default: ;
    }
#endif
#endif
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
    VertexInput input (i1, w, v, vpp, i_in, spin, channel);
    gammaRb = vertex[0].gammaRb(input);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex[0].avertex.valsmooth(k2, input, vertex[0].tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].avertex.valsmooth(k3, input, vertex[0].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex[0].pvertex.valsmooth(k2, input);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].pvertex.valsmooth(k3, input);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex[0].tvertex.valsmooth(k2, input, vertex[0].avertex);
#endif
#if DIAG_CLASS >=3
            K3 = vertex[0].tvertex.valsmooth(k3, input, vertex[0].avertex);
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
    VertexInput input (i3, w, vpp, vp, i_in, spin, channel);
    gammaRb = vertex[0].gammaRb(input);
#endif

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex[0].avertex.valsmooth(k2b, input, vertex[0].tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].avertex.valsmooth(k3, input, vertex[0].tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex[0].pvertex.valsmooth(k2b, input);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].pvertex.valsmooth(k3, input);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex[0].tvertex.valsmooth(k2b, input, vertex[0].avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = vertex[0].tvertex.valsmooth(k3, input, vertex[0].avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}



#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
