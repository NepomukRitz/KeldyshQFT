//
// Created by E.Walter on 7/31/19.
//

#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

#include "frequency_grid.h"

//TODO: naming??

/******************CLASS FOR SELF ENERGY *************/
template <typename Q>
class SelfEnergy{
    vec<Q> Sigma =  vec<Q> (2*nSE); // factor 2 for Keldysh components: Sigma^R, Sigma^K TODO: check
public:
    void setself(int, int, Q);
    Q sval(int, int);
    Q svalsmooth(int, double);
    friend SelfEnergy operator+(const SelfEnergy& self1, const SelfEnergy& self2);
    friend SelfEnergy operator+=(const SelfEnergy& self1,const SelfEnergy& self2);
    friend SelfEnergy operator*(Q alpha, const SelfEnergy& self1);
    friend SelfEnergy operator*(const SelfEnergy& self1, Q alpha);
};



/*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
template <typename Q>
Q SelfEnergy<Q>::sval(int iK, int i){
    return Sigma[iK*nSE + i];
}

template <typename Q>
Q SelfEnergy<Q>::svalsmooth(int iK, double w){//smoothly interpolates for values between discrete frequency values of mesh

    int W = fconv(w);

    double x1 = ffreqs[W];
    double x2 = ffreqs[W+1];
    double xd = (w-x1)/(x2-x1);

    Q f1 = sval(iK, W);
    Q f2 = sval(iK, W+1);

    return (1.-xd)*f1 + xd*f2;
}

template <typename Q>
void SelfEnergy<Q>::setself(int iK, int i, Q val){
    Sigma[iK*nSE + i] = val;
}

//operators for self energy
template <typename Q>
SelfEnergy<Q> operator*(Q alpha, const SelfEnergy<Q>& self1){//product operator overloading
    SelfEnergy<Q> self2;
    self2.Sigma = self1.Sigma * alpha;
    return self2;
}
template <typename Q>
SelfEnergy<Q> operator*(const SelfEnergy<Q>& self1, Q alpha){//product operator overloading
    SelfEnergy<Q> self2;
    self2.Sigma = self1.Sigma * alpha;
    return self2;
}
template <typename Q>
SelfEnergy<Q> operator+(const SelfEnergy<Q>& self1, const SelfEnergy<Q>& self2){//sum operator overloading
    SelfEnergy<Q> self3;
    self3.Sigma = self1.Sigma + self2.Sigma;
    return self3;
}
template <typename Q>
SelfEnergy<Q> operator+=(const SelfEnergy<Q>& self1, const SelfEnergy<Q>& self2){//sum operator overloading
    SelfEnergy<Q> self3;
    self3.Sigma = self1.Sigma + self2.Sigma;
    return self3;
}






#endif //KELDYSH_MFRG_SELFENERGY_H
