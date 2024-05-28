//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP

#include "../parameters/master_parameters.hpp"
#include "../utilities/math_utils.hpp"
#include "Keldysh_symmetries.hpp"
#include "../data_structures.hpp"


void switch_channel(IndicesSymmetryTransformations& indices);

/// exchange incoming legs
template <bool is_SBE_lambda> void T1 (IndicesSymmetryTransformations& indices) {
    if constexpr (KELDYSH_FORMALISM) {
        std::array<int,4> alphas_tmp = iK_to_alphas(indices.iK);
        std::swap(alphas_tmp[2],alphas_tmp[3]);
        indices.iK = alphas_to_iK(alphas_tmp);
    }
    if constexpr(!is_SBE_lambda)
    {
        indices.prefactor *= -1.;
    }
    //else {
    //    if (indices.channel_rvert == 'p') indices.prefactor *= -1.;
    //}
    indices.conjugate ^= false;

    if(indices.channel_rvert == 'p'){
        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= 1.;
            indices.v2 *= -1.;
            if (!KELDYSH && !ZERO_T){
                freqType rounding_correction = signFlipCorrection_MF(indices.w);  // correction due to rounding towards Matsubara frequencies
                indices.v2 += rounding_correction;
                // correction due to rounding towards Matsubara frequencies
            }
        }
    }
    else {
        indices.asymmetry_transform ^= true;

        indices.w  *= -1.;
        if (MAX_DIAG_CLASS > 1) {
            freqType temp = indices.v1;
            indices.v1 = indices.v2;
            indices.v2 = temp;
        }
    }


    switch_channel(indices);
    indices.spin = 1 - indices.spin;

}
/// exchange outgoing lines
template <bool is_SBE_lambda> void T2 (IndicesSymmetryTransformations& indices) {
    if constexpr (KELDYSH_FORMALISM) {
        std::array<int,4> alphas_tmp = iK_to_alphas(indices.iK);
        std::swap(alphas_tmp[0],alphas_tmp[1]);
        indices.iK = alphas_to_iK(alphas_tmp);
    }
    if constexpr(!is_SBE_lambda)
    {
        indices.prefactor *= -1.;
    }
    //else {
    //    if (indices.channel_rvert == 'p') indices.prefactor *= -1.;
    //}
    indices.conjugate ^= false;

    if(indices.channel_rvert == 'p'){
        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1.;
            indices.v2 *= 1.;
            if (!KELDYSH && !ZERO_T) {
                double rounding_correction = signFlipCorrection_MF(indices.w);  // correction due to rounding towards Matsubara frequencies
                indices.v1 += rounding_correction;
                // correction due to rounding towards Matsubara frequencies
            }
        }
    }
    switch_channel(indices);
    indices.spin = 1 - indices.spin;
}
template <bool is_SBE_lambda> void T3 (IndicesSymmetryTransformations& indices) {
    T1<is_SBE_lambda>(indices);
    T2<is_SBE_lambda>(indices);
}
template <bool is_SBE_lambda> void TC (IndicesSymmetryTransformations& indices) {
    if constexpr (KELDYSH_FORMALISM) {
        std::array<int,4> alphas_tmp = iK_to_alphas(indices.iK);
        std::swap(alphas_tmp[0],alphas_tmp[2]);
        std::swap(alphas_tmp[1],alphas_tmp[3]);
        indices.iK = alphas_to_iK(alphas_tmp);
    }
    indices.conjugate ^= true;


    if (indices.channel_rvert == 't'){ //TC acts differently on t, not on p!!
        indices.w  *= -1.;
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


    if constexpr(KELDYSH and CONTOUR_BASIS != 1)
    {

        if constexpr(is_SBE_lambda) {
            if ( isInList(indices.iK, odd_Keldysh)) indices.prefactor *= -1.;
        }
        else  {
            if (!isInList(indices.iK, odd_Keldysh)) indices.prefactor *= -1.;
        }
    }
    else if constexpr(KELDYSH and CONTOUR_BASIS == 1)
    {

        if constexpr(!is_SBE_lambda) {indices.prefactor *= -1;}
    }
    else{
        indices.w *= -1;
        indices.v1 *= -1;
        indices.v2 *= -1;
        indices.iw_r = nBOS3 - 1 - indices.iw_r; // corresponds to "w *= -1"
        if (!ZERO_T){// Matsubara T>0
            double rounding_correction = signFlipCorrection_MF(indices.w);  // correction due to rounding towards Matsubara frequencies
            indices.v1 += rounding_correction;
            indices.v2 += rounding_correction;
        }
    }
}
template <bool is_SBE_lambda> void Tph (IndicesSymmetryTransformations& indices) {
    if (KELDYSH && PARTICLE_HOLE_SYMMETRY){ // Used only for Keldysh calculations with particle-hole symmetry
        indices.conjugate ^= true;

            if constexpr(CONTOUR_BASIS != 1)
            {
                if constexpr(is_SBE_lambda) {
                    if ( isInList(indices.iK, odd_Keldysh)) indices.prefactor *= -1.;
                }
                else  {
                    if (!isInList(indices.iK, odd_Keldysh)) indices.prefactor *= -1.;
                }
            }
            else {
                if constexpr(!is_SBE_lambda) {indices.prefactor *= -1.;}
            }


        indices.w *= -1;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1;
            indices.v2 *= -1;
        }
    }
}
template <bool is_SBE_lambda> void TR (IndicesSymmetryTransformations& indices) {
    if (!KELDYSH){ // used only for Matsubara calculations
        indices.conjugate ^= true;
        indices.w *= -1;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= -1;
            indices.v2 *= -1;
            if (!ZERO_T) { // Matsubara T>0
                double rounding_correction = signFlipCorrection_MF(indices.w);
                indices.v1 += rounding_correction;
                indices.v2 += rounding_correction;
            }
        }
    }
}
template <bool is_SBE_lambda> void Ti (IndicesSymmetryTransformations& indices, int i)  {
    if (i == 0) return;
    switch (i) {
        case 1:
            T1<is_SBE_lambda>(indices);
            break;
        case 2:
            T2<is_SBE_lambda>(indices);
            break;
        case 3:
            T3<is_SBE_lambda>(indices);
            break;
        case 4:
            TC<is_SBE_lambda>(indices);
            break;
        case 14:
            T1<is_SBE_lambda>(indices);
            TC<is_SBE_lambda>(indices);
            break;
        case 24:
            T2<is_SBE_lambda>(indices);
            TC<is_SBE_lambda>(indices);
            break;
        case 41:
            TC<is_SBE_lambda>(indices);
            T1<is_SBE_lambda>(indices);
            break;
        case 43:
            TC<is_SBE_lambda>(indices);
            T3<is_SBE_lambda>(indices);
            break;
        case 34:
            TC<is_SBE_lambda>(indices);
            T3<is_SBE_lambda>(indices);
            break;
#if KELDYSH_FORMALISM and PARTICLE_HOLE_SYMM
                case 6:
                    Tph<is_SBE_lambda>(indices);
                break;
                case 16:
                    T1<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 26:
                    T2<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 36:
                    T3<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 46:
                    TC<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 146:
                    T1<is_SBE_lambda>(indices);
                TC<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 246:
                    T2<is_SBE_lambda>(indices);
                TC<is_SBE_lambda>(indices);
                Tph<is_SBE_lambda>(indices);
                break;
                case 346:
                    Tph<is_SBE_lambda>(indices);
                TC<is_SBE_lambda>(indices);
                T3<is_SBE_lambda>(indices);
                break;
#endif
            if (!KELDYSH) {
                case 7:
                    TR<is_SBE_lambda>(indices);
                break;
                case 37:
                    TR<is_SBE_lambda>(indices);
                T3<is_SBE_lambda>(indices);
                break;
                case 47:
                    TC<is_SBE_lambda>(indices);
                TR<is_SBE_lambda>(indices);
                break;
                case 347:
                    TR<is_SBE_lambda>(indices);
                TC<is_SBE_lambda>(indices);
                T3<is_SBE_lambda>(indices);
                break;
            }
        default:
            utils::print("The Transformation ", i, " in the symmetry table is not covered in Ti! Abort."); assert(false);
    }
}

// test cases:
// 1.) obvious checks for some concrete values
// 2.) group structure of transformations: T1TC == TCT2, T1T2 == T2T1, T1TC != TCT1, etc., for all channels


#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP