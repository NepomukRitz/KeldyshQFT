//
// Created by Sa.Aguirre on 8/16/19.
//

#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_H

#include <array>

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
int T_3_Keldysh(int iK)
{
    int alpha1p, alpha2p, alpha1, alpha2, iKp;
    alpha2 = iK % 2;
    alpha1 = (iK % 4) / 2;
    alpha2p = (iK % 8) / 4;
    alpha1p = (iK % 16) / 8;

    iKp = alpha2p*8 + alpha1p*4 + alpha2*2 + alpha1;

    return iKp;
}

/*T_C interchanges incoming with outgoing legs i.e. complex conjugation*/
int T_C_Keldysh(int iK)
{
    int alpha1p, alpha2p, alpha1, alpha2, iKp;
    alpha2 = iK % 2;
    alpha1 = (iK % 4) / 2;
    alpha2p = (iK % 8) / 4;
    alpha1p = (iK % 16) / 8;

    iKp = alpha1*8 + alpha2*4 + alpha1p*2 + alpha2p;

    return iKp;
}

bool isInList (int iK, vector<int> list)
{
    bool resp = false;
    int i=0;
    while(!resp && i<list.size())
    {
        if(iK == list[i])
        {
            resp = true;
        }
        i++;
    }
    return resp;
}


//*The next functions are the group of actions that transform related Keldysh components within the vertex into one another.
// * In lack of a better nomenclature, name them A, since they switch alphas...*/
//int A1_K1a(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p  )%2;
//    index1 = (alpha1   )%2;
//    index0 = (alpha2+1 )%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A2_K1a(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p  )%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2   )%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A3_K1a(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2 +1)%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//
//
//int A1_K1p(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1   )%2;
//    index0 = (alpha2   )%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A2_K1p(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p  )%2;
//    index2 = (alpha2p  )%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2 +1)%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A3_K1p(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2 +1)%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//
//
//
//
//int A1_K1t(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p  )%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1   )%2;
//    index0 = (alpha2 +1)%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A2_K1t(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p  )%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2   )%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}
//
//int A3_K1t(int iK)
//{
//    int alpha1p,alpha2p, alpha1, alpha2;
//    int index0, index1, index2, index3;
//
//    alpha2 = iK%2;
//    alpha1 = (iK%4)/2;
//    alpha2p = (iK%8)/4;
//    alpha1p = (iK%16)/8;
//
//    index3 = (alpha1p+1)%2;
//    index2 = (alpha2p+1)%2;
//    index1 = (alpha1 +1)%2;
//    index0 = (alpha2 +1)%2;
//
//    return (index3)*8 + (index2)*4 + (index1)*2 + (index0);
//}

/* This function converts indices in the range 0...5 to the actual Keldysh index they correspond to
 * Rule: {0,1,2,3,4,5} -> {0,1,3,5,6,7}
 * The components 0,1,3,5,6 and 7 are the chosen reference components, numerated in the 0...15 "convention"*/
int convertToRealIndex(int index)
{
    if(index ==0 || index == 1)
        return index;
    else if(index == 2)
        return 3;
    else
        return index + 2;
}

int convertToIndepIndex(int iK)
{
    if(iK ==0 || iK ==1)
        return iK;
    else if(iK==3)
        return 2;
    else if(iK==5)
        return 3;
    else if(iK ==7)
        return 4;
}

/*This function returns the values of the 4 alphas for a given index in the 0...15 set */
tuple<int, int, int, int> alphas(int index)
{
    int alpha1p, alpha2p, alpha1, alpha2;

    alpha2 = index % 2+1;
    alpha1 = (index % 4) / 2+1;
    alpha2p = (index % 8) / 4+1;
    alpha1p = (index % 16) / 8+1;

    return make_tuple(alpha1p, alpha2p, alpha1, alpha2);
}





#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
