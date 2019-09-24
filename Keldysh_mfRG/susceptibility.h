//
// Created by Sa.Aguirre on 8/20/19.
//

#ifndef KELDYSH_MFRG_SUSCEPTIBILITY_H
#define KELDYSH_MFRG_SUSCEPTIBILITY_H

#include "data_structures.h"
#include "parameters.h"
#include "frequency_grid.h"

#ifdef SUSC

/******************CLASS FOR SUSCEPTIBILITY *************/
template <typename Q>
class Susc{
public:
    vec<Q> Susceptibility =  vec<Q> (2*nSUSC); // factor 2 for Keldysh components: Susc^R, Susc^K. Susc^A = [Sus^R] dagger

    void setsus(int, int, Q);
    Q susval(int, int);
    Q susvalsmooth(int, double);
//    friend Susc operator+(const Susc& sus1, const Susc& sus2);
//    friend Susc operator+=(const Susc& sus1,const Susc& sus2);
//    friend Susc operator*(Q alpha, const Susc& sus1);
//    friend Susc operator*(const Susc& sus1, Q alpha);
};



/*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
template <typename Q>Q Susc<Q>::susval(int iK, int i){
    return Susceptibility[iK*nSE + i];
}
template <typename Q>Q Susc<Q>::susvalsmooth(int iK, double w){//smoothly interpolates for values between discrete frequency values of mesh

    int W = fconv(w);

    double x1 = ffreqs[W];
    double x2 = ffreqs[W+1];
    double xd = (w-x1)/(x2-x1);

    Q f1 = susval(iK, W);
    Q f2 = susval(iK, W+1);

    return (1.-xd)*f1 + xd*f2;
}
template <typename Q>void Susc<Q>::setsus(int iK, int i, Q val){
    Susceptibility[iK*nSE + i] = val;
}

//operators for susceptibility
template <typename Q>Susc<Q> operator+(const Susc<Q>& sus1, const Susc<Q>& sus2){//sum operator overloading
    Susc<Q> sus3;
    sus3.Susceptibility = sus1.Susceptibility + sus2.Susceptibility;
    return sus3;
}
template <typename Q>Susc<Q> operator+=(const Susc<Q>& sus1, const Susc<Q>& sus2){//sum operator overloading
    Susc<Q> sus3;
    sus3.Susceptibility = sus1.Susceptibility + sus2.Susceptibility;
    return sus3;
}
template <typename Q>Susc<Q> operator*(Q alpha, const Susc<Q>& sus1){//product operator overloading
    Susc<Q> sus2;
    sus2.Susceptibility = sus1.Susceptibility * alpha;
    return sus2;
}
template <typename Q>Susc<Q> operator*(const Susc<Q>& sus1, Q alpha){//product operator overloading
    Susc<Q> sus2;
    sus2.Susceptibility = sus1.Susceptibility * alpha;
    return sus2;
}
template <typename Q>Susc<Q> operator*(double alpha, const Susc<Q>& sus1){//product operator overloading
    Susc<Q> sus2;
    sus2.Susceptibility = sus1.Susceptibility * alpha;
    return sus2;
}
template <typename Q>Susc<Q> operator*(const Susc<Q>& sus1, double alpha){//product operator overloading
    Susc<Q> sus2;
    sus2.Susceptibility = sus1.Susceptibility * alpha;
    return sus2;
}
#endif


#endif //KELDYSH_MFRG_SUSCEPTIBILITY_H
