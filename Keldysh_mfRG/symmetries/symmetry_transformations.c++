#include "symmetry_transformations.hpp"

void switch_channel(IndicesSymmetryTransformations& indices) {
    switch (indices.channel) {
        case 'a':
            indices.channel = 't';
            break;
        case 't':
            indices.channel = 'a';
            break;
        default:;
    }
    if (INTERPOL2D_FOR_K3) {
        switch (indices.channel_bubble) {
            case 'a':
                indices.channel_bubble = 't';
                break;
            case 't':
                indices.channel_bubble = 'a';
                break;
            default: ;
        }
    }
}

void T1 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
        indices.w  *= 1.;
        if (MAX_DIAG_CLASS > 1) {
            indices.v1 *= 1.;
            indices.v2 *= -1.;
            if (!KELDYSH && !ZERO_T){
                double rounding_correction = signFlipCorrection_MF(indices.w);  // correction due to rounding towards Matsubara frequencies
                indices.v2 += rounding_correction;
                // correction due to rounding towards Matsubara frequencies
            }
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

    if (INTERPOL2D_FOR_K3) {
        if(indices.channel_bubble != 'p'){
            indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
        }
    }

    switch_channel(indices);
    indices.spin = 1 - indices.spin;

}

void T2 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
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
    else {
//        indices.w  *= 1.;
//        indices.v1 *= 1.;
//        indices.v2 *= 1.;
    }
    switch_channel(indices);
    indices.spin = 1 - indices.spin;
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

    if (INTERPOL2D_FOR_K3) {
        if (indices.channel_bubble == 't') {
            indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
        }
    }

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
            double rounding_correction = signFlipCorrection_MF(indices.w);  // correction due to rounding towards Matsubara frequencies
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

        if (INTERPOL2D_FOR_K3) {
            indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
        }
    }
}

void TR (IndicesSymmetryTransformations& indices){
    if (!KELDYSH){ // used only for Matsubara calculations
        indices.conjugate ^= true;
        indices.w *= -1;
        if (INTERPOL2D_FOR_K3) {
            indices.iw = nBOS3 - 1 - indices.iw; // corresponds to "w *= -1"
        }
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
