//
// Created by E.Walter on 7/31/19.
//

#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H


#include "parameters.h"

//TODO: naming??

/******************CLASS FOR SELF ENERGY *************/
template <typename Q>
class self{
    vec<Q> selfenergy =  vec<Q> (2*nSE); // factor 2 for Keldysh components: Sigma^R, Sigma^K TODO: check
public:
    void setself(int, int, Q);
    Q sval(int, int);
    Q svalsmooth(int, double);
    friend self operator+(const self& self1, const self& self2);
    friend self operator+=(const self& self1,const self& self2);
    friend self operator*(Q alpha, const self& self1);
    friend self operator*(const self& self1, Q alpha);
};



/*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
template <typename Q>
Q self<Q>::sval(int iK, int i){
    return selfenergy[iK*nSE + i];
}

template <typename Q>
Q self<Q>::svalsmooth(int iK, double w){//smoothly interpolates for values between discrete frequency values of mesh
    Q value;
    int W = fconv(w); // TODO: define fconv
    value += ((selfenergy[iK*nSE+W]*(ffreqs[W+1]-w)+selfenergy[iK*nSE+W+1]*(-ffreqs[W]+w))/(ffreqs[W+1]-ffreqs[W])); // TODO: make it similar to the version for vertex
    return value;
}

template <typename Q>
void self<Q>::setself(int iK, int i, Q val){
    selfenergy[iK*nSE + i] = val;
}

//operators for self energy
self operator*(double alpha, const self& self1){//product operator overloading
    self self2;
    self2.selfenergy = self1.selfenergy * alpha;
    return self2;
}
self operator*(const self& self1, double alpha){//product operator overloading
    self self2;
    self2.selfenergy = self1.selfenergy * alpha;
    return self2;
}
self operator+(const self& self1, const self& self2){//sum operator overloading
    self self3;
    self3.selfenergy = self1.selfenergy + self2.selfenergy;
    return self3;
}
self operator+=(const self& self1, const self& self2){//sum operator overloading
    self self3;
    self3.selfenergy = self1.selfenergy + self2.selfenergy;
    return self3;
}






#endif //KELDYSH_MFRG_SELFENERGY_H
