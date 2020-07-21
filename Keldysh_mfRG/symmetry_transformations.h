//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H

#include "parameters.h"

struct IndicesSymmetryTransformations{
    int iK;
    double prefactor;
    bool conjugate;
    bool transform;
    double w, v1, v2; int i_in;

    IndicesSymmetryTransformations(int iK_in, double w_in, double v1_in, double v2_in, int i_in_in)
            : iK(iK_in), w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in)
    {
        prefactor = 1.;
        conjugate = false;
        transform = false;
    }
};

void T1 (IndicesSymmetryTransformations& indicesSymmetryTransformations, const char channel){

    indicesSymmetryTransformations.prefactor *= -1.;
    indicesSymmetryTransformations.conjugate ^= false;
    indicesSymmetryTransformations.transform ^= true;

    if(channel == 'p'){
        indicesSymmetryTransformations.w  *= 1.;
#if DIAG_CLASS>1
        indicesSymmetryTransformations.v1 *= 1.;
        indicesSymmetryTransformations.v2 *= -1.;
#endif
    }
    else {
        indicesSymmetryTransformations.w  *= -1.;
#if DIAG_CLASS>1
        double temp = indicesSymmetryTransformations.v1;
        indicesSymmetryTransformations.v1 = indicesSymmetryTransformations.v2;
        indicesSymmetryTransformations.v2 = temp;
#endif
    }
}

void T2 (IndicesSymmetryTransformations& indicesSymmetryTransformations, const char channel){

    indicesSymmetryTransformations.prefactor *= -1.;
    indicesSymmetryTransformations.conjugate ^= false;
    indicesSymmetryTransformations.transform ^= true;

    if(channel == 'p'){
        indicesSymmetryTransformations.w  *= 1.;
#if DIAG_CLASS>1
        indicesSymmetryTransformations.v1 *= -1.;
        indicesSymmetryTransformations.v2 *= 1.;
#endif
    }
    else {
//        indicesSymmetryTransformations.w  *= 1.;
//        indicesSymmetryTransformations.v1 *= 1.;
//        indicesSymmetryTransformations.v2 *= 1.;
    }
}

void T3 (IndicesSymmetryTransformations& indicesSymmetryTransformations, const char channel){
    T1(indicesSymmetryTransformations, channel);
    T2(indicesSymmetryTransformations, channel);
}

void TC (IndicesSymmetryTransformations& indicesSymmetryTransformations, const char channel){

    if(isInList(indicesSymmetryTransformations.iK, odd_Keldysh))
        indicesSymmetryTransformations.prefactor *= 1.;
    else
        indicesSymmetryTransformations.prefactor *= -1.;

    indicesSymmetryTransformations.conjugate ^= true;
    indicesSymmetryTransformations.transform ^= false;

    if (channel == 't'){ //TC acts differently on t, not on p!!
        indicesSymmetryTransformations.w  *= -1.;
//        indicesSymmetryTransformations.v1 *= 1.;
//        indicesSymmetryTransformations.v2 *= 1.;
    }
    else{
        indicesSymmetryTransformations.w  *= 1.;
#if DIAG_CLASS>1
        double temp = indicesSymmetryTransformations.v1;
        indicesSymmetryTransformations.v1 = indicesSymmetryTransformations.v2;
        indicesSymmetryTransformations.v2 = temp;
#endif
    }
}

#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_H