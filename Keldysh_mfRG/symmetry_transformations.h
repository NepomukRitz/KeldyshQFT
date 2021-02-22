//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H

#include "parameters.h"
#include "Keldysh_symmetries.h"

/**
 * specifies the diagrammatic contribution ( + modifications by prefactor, complex conjugation)
 */
struct IndicesSymmetryTransformations{
    int iK;
    double prefactor = 1.;
    bool conjugate = false;
    bool asymmetry_transform = false;
    double w, v1, v2; int i_in;
    char channel;

    IndicesSymmetryTransformations(int iK_in, double w_in, double v1_in, double v2_in, int i_in_in, char channel_in)
            : iK(iK_in), w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in), channel(channel_in)
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
}

void T1 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
        indices.w  *= 1.;
#if DIAG_CLASS>1
        indices.v1 *= 1.;
        indices.v2 *= -1.;
#endif
    }
    else {
        indices.asymmetry_transform ^= true;

        indices.w  *= -1.;
#if DIAG_CLASS>1
        double temp = indices.v1;
        indices.v1 = indices.v2;
        indices.v2 = temp;
#endif
    }
    switch_channel(indices);
}

void T2 (IndicesSymmetryTransformations& indices){

    indices.prefactor *= -1.;
    indices.conjugate ^= false;

    if(indices.channel == 'p'){
        indices.w  *= 1.;
#if DIAG_CLASS>1
        indices.v1 *= -1.;
        indices.v2 *= 1.;
#endif
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
#if DIAG_CLASS>1
        double temp = indices.v1;
        indices.v1 = indices.v2;
        indices.v2 = temp;
#endif
    }

#ifdef KELDYSH_FORMALISM
    if(isInList(indices.iK, odd_Keldysh))
        indices.prefactor *= 1.;
    else
        indices.prefactor *= -1.;
#else
    indices.w *= -1;
    indices.v1 *= -1;
    indices.v2 *= -1;
#endif
}

#ifndef KELDYSH_FORMALISM
void TR (IndicesSymmetryTransformations& indices){
    indices.conjugate ^= true;
    indices.w *= -1;
#if DIAG_CLASS > 1
    indices.v1 *= -1;
    indices.v2 *= -1;
#endif
}
#endif

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
#ifndef KELDYSH_FORMALISM
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
#endif
        default:
            cout << "A Transformation in the symmetry table is not covered in Ti!"
            ;
    }
}

// test cases:
// 1.) obvious checks for some concrete values
// 2.) group structure of transformations: T1TC == TCT2, T1T2 == T2T1, T1TC != TCT1, etc., for all channels



#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H