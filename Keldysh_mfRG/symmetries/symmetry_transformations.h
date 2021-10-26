//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H

#include "../parameters/master_parameters.h"
#include "../utilities/math_utils.h"
#include "Keldysh_symmetries.h"

/**
 * specifies the diagrammatic contribution ( + modifications by prefactor, complex conjugation)
 */
struct IndicesSymmetryTransformations{
    int iK;
    double prefactor = 1.; // fermionic sign factor; comes in effect for T1, T2 (and sometimes Tc)
    bool conjugate = false;
    bool asymmetry_transform = false;
    int iw;
    double w, v1, v2; int i_in;
    K_class kClass_aim;     // we only distinguish kClass==k3 from kClass!=k3 --> important for interpolation in K3
    char channel;
    char channel_bubble;
    char channel_parametrization = channel; // W.r.t. which channel is the vertex parametrized? Used for the Hubbard model.

    IndicesSymmetryTransformations(int iK_in, double w_in, double v1_in, double v2_in, int i_in_in, char channel_in, K_class k_in, int iw_in, char channel_bubble_in)
            : iK(iK_in), w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in), channel(channel_in), kClass_aim(k_in), iw(iw_in), channel_bubble(channel_bubble_in)
    {}

    IndicesSymmetryTransformations(VertexInput input, char channel_in)
            : iK(input.iK), w(input.w), v1(input.v1), v2(input.v2), i_in(input.i_in), channel(channel_in), kClass_aim(input.kClass_aim), iw(input.iw_r), channel_bubble(input.channel)
    {}
};

void switch_channel(IndicesSymmetryTransformations& indices) {
    switch (indices.channel) {
        case 'a':
            indices.channel = 't';
            break;
        case 't':
            indices.channel = 'a';
            break;
        default: ;
    }
#ifdef INTERPOL2D_FOR_K3
    switch (indices.channel_bubble) {
        case 'a':
            indices.channel_bubble = 't';
            break;
        case 't':
            indices.channel_bubble = 'a';
            break;
        default: ;
    }
#endif
}

void T1 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= 1.;
            indices.v2 *= -1.;
            if (!KELDYSH && !ZERO_T)
                indices.v2 += floor2bfreq(indices.w / 2) - ceil2bfreq(indices.w / 2);
            // correction due to rounding towards Matsubara frequencies
        }
    }
    else {
        indices.asymmetry_transform ^= true;

        indices.w  *= -1.;
        if (MAX_DIAG_CLASS > 1) {
            double temp = indices.v1;
            indices.v1 = indices.v2;
            indices.v2 = temp;
        }
    }

#ifdef INTERPOL2D_FOR_K3
    if(indices.channel_bubble != 'p'){
        indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
    }
#endif
    switch_channel(indices);
}

void T2 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1.;
            indices.v2 *= 1.;
            if (!KELDYSH && !ZERO_T)
                indices.v1 += floor2bfreq(indices.w / 2) - ceil2bfreq(indices.w / 2);
            // correction due to rounding towards Matsubara frequencies
        }
    }
    else {
//        indices.w  *= 1.;
//        indices.v1 *= 1.;
//        indices.v2 *= 1.;
    }
    switch_channel(indices);
}

void T3 (IndicesSymmetryTransformations& indices){
    T1(indices);
    T2(indices);
}

void TC (IndicesSymmetryTransformations& indices){

    indices.conjugate ^= true;

    if (indices.channel == 't'){ //TC acts differently on t, not on p!!
        indices.w  *= -1.;
//        indices.v1 *= 1.;
//        indices.v2 *= 1.;
    }
    else{
        indices.asymmetry_transform ^= true;

        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            double temp = indices.v1;
            indices.v1 = indices.v2;
            indices.v2 = temp;
        }
    }

#ifdef INTERPOL2D_FOR_K3
    if (indices.channel_bubble == 't') {
        indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
    }
#endif

    if (KELDYSH){
        if(isInList(indices.iK, odd_Keldysh))
            indices.prefactor *= 1.;
        else
            indices.prefactor *= -1.;
    }
    else{
        indices.w *= -1;
        indices.v1 *= -1;
        indices.v2 *= -1;
        indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
        if (!ZERO_T){// Matsubara T>0
            double rounding_correction = floor2bfreq(indices.w/2) -  ceil2bfreq(indices.w/2);  // correction due to rounding towards Matsubara frequencies
            indices.v1 += rounding_correction;
            indices.v2 += rounding_correction;
        }
    }
}

void Tph (IndicesSymmetryTransformations& indices){
    if (KELDYSH && PARTICLE_HOLE_SYMMETRY){ // Used only for Keldysh calculations with particle-hole symmetry
        indices.conjugate ^= true;
        if (isInList(indices.iK, odd_Keldysh))
            indices.prefactor *= 1.;
        else
            indices.prefactor *= -1.;

        indices.w *= -1;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1;
            indices.v2 *= -1;
        }

#ifdef INTERPOL2D_FOR_K3
        indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
#endif
    }
}

void TR (IndicesSymmetryTransformations& indices){
    if (!KELDYSH){ // used only for Matsubara calculations
        indices.conjugate ^= true;
        indices.w *= -1;
#ifdef INTERPOL2D_FOR_K3
        indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
#endif
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1;
            indices.v2 *= -1;
            if (!ZERO_T) { // Matsubara T>0
                double rounding_correction = floor2bfreq(indices.w / 2) - ceil2bfreq(
                        indices.w / 2);  // correction due to rounding towards Matsubara frequencies
                indices.v1 += rounding_correction;
                indices.v2 += rounding_correction;
            }
        }
    }
}

/**
 * Define symmetry transformations.
 * The function modifies the indices according to the transformation T_i.
 * @param indices     : specifies the diagrammatic contribution ( + modifications by prefactor, complex conjugation)
 * @param i           : specities the transformation T_i
 */
void Ti (IndicesSymmetryTransformations& indices, const int i) {
    if (i == 0) return;
    switch (i) {
        case 1:
            T1(indices);
            break;
        case 2:
            T2(indices);
            break;
        case 3:
            T3(indices);
            break;
        case 4:
            TC(indices);
            break;
        case 14:
            T1(indices);
            TC(indices);
            break;
        case 41:
            TC(indices);
            T1(indices);
            break;
        case 43:
            TC(indices);
            T3(indices);
            break;
        case 34:
            TC(indices);
            T3(indices);
            break;
    if (KELDYSH and PARTICLE_HOLE_SYMMETRY) {
                case 6:
                    Tph(indices);
            break;
        case 36:
            Tph(indices);
            T3(indices);
            break;
        case 46:
            Tph(indices);
            TC(indices);
            break;
        case 346:
            Tph(indices);
            TC(indices);
            T3(indices);
            break;
    }
    if (!KELDYSH) {
        case 7:
            TR(indices);
        break;
        case 37:
            TR(indices);
        T3(indices);
        break;
        case 47:
            TC(indices);
        TR(indices);
        break;
        case 347:
            TR(indices);
        TC(indices);
        T3(indices);
        break;
    }
        default:
            print("A Transformation in the symmetry table is not covered in Ti! Abort."); assert(false);
    }
}

// test cases:
// 1.) obvious checks for some concrete values
// 2.) group structure of transformations: T1TC == TCT2, T1T2 == T2T1, T1TC != TCT1, etc., for all channels



#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H