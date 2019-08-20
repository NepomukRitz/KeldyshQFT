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


/*The next functions are the group of actions that transform related Keldysh components within the vertex into one another.
 * In lack of a better nomenclature, name them A, since they switch alphas...*/
int A1(int iK)
{
    int alpha1p,alpha2p, alpha1, alpha2;
    int index0, index1, index2, index3;

    alpha2 = iK%2-1;
    alpha1 = (iK%4)/2-1;
    alpha2p = (iK%8)/4-1;
    alpha1p = (iK%16)/8-1;

    index3 = (alpha1p+1)%2;
    index2 = (alpha2p  )%2;
    index1 = (alpha1   )%2;
    index0 = (alpha2+1) %2;

    return (index3+1)*8 + (index2+1)*4 + (index1+1)*2 + (index0+1);
}

int A2(int iK)
{
    int alpha1p,alpha2p, alpha1, alpha2;
    int index0, index1, index2, index3;

    alpha2 = iK%2-1;
    alpha1 = (iK%4)/2-1;
    alpha2p = (iK%8)/4-1;
    alpha1p = (iK%16)/8-1;

    index3 = (alpha1p  )%2;
    index2 = (alpha2p+1)%2;
    index1 = (alpha1 +1)%2;
    index0 = (alpha2   )%2;

    return (index3+1)*8 + (index2+1)*4 + (index1+1)*2 + (index0+1);
}

int A3(int iK)
{
    int alpha1p,alpha2p, alpha1, alpha2;
    int index0, index1, index2, index3;

    alpha2 = iK%2-1;
    alpha1 = (iK%4)/2-1;
    alpha2p = (iK%8)/4-1;
    alpha1p = (iK%16)/8-1;

    index3 = (alpha1p+1)%2;
    index2 = (alpha2p+1)%2;
    index1 = (alpha1 +1)%2;
    index0 = (alpha2 +1)%2;

    return (index3+1)*8 + (index2+1)*4 + (index1+1)*2 + (index0+1);
}


#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
