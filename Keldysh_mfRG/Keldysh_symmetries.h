//
// Created by Sa.Aguirre on 8/16/19.
//

#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_H

/*T_1 switches the incoming legs*/
int T_1_Keldysh(int iK)
{
    int alpha1p, alpha2p, alpha1, alpha2, iKp;

    alpha2 = iK%2;
    alpha1 = (iK%4)/2;
    alpha2p = (iK%8)/4;
    alpha1p = (iK%16)/8;

    iKp = alpha1p*8 + alpha2p*4 + alpha2*2 + alpha1;

    return iKp;
}

/* T2 switches the outgoing legs */
int T_2_Keldysh(int iK)
{
    int alpha1p, alpha2p, alpha1, alpha2, iKp;

    alpha2 = iK%2;
    alpha1 = (iK%4)/2;
    alpha2p = (iK%8)/4;
    alpha1p = (iK%16)/8;

    iKp = alpha2p*8 + alpha1p*4 + alpha1*2 + alpha2;

    return iKp;
}

/* T3 switches both incoming and outgoing legs */
int T_3_Keldysh(int iK) {
    int alpha1p, alpha2p, alpha1, alpha2, iKp;

    alpha2 = iK % 2;
    alpha1 = (iK % 4) / 2;
    alpha2p = (iK % 8) / 4;
    alpha1p = (iK % 16) / 8;

    iKp = alpha2p * 8 + alpha1p * 4 + alpha2 * 2 + alpha1;

    return iKp;
}

#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
