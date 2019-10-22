//
// Created by E.Walter on 7/31/19.
//

#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H


/******************CLASS FOR SELF ENERGY *************/
template <typename Q>
class SelfEnergy{
    vec<Q> Sigma =  vec<Q> (2*nSE); // factor 2 for Keldysh components: Sigma^R, Sigma^K
public:
    Q sval(int, int);
    Q svalsmooth(int, double);
    void setself(int, int, Q);

//operators for self energy

    SelfEnergy<Q> operator+(const SelfEnergy<Q>& self1){//sum operator overloading
        this->Sigma + self1.Sigma;
        return *this;
    }
    SelfEnergy<Q> operator+=(const SelfEnergy<Q>& self1){//sum operator overloading
        this->Sigma += self1.Sigma;
        return *this;
    }
    SelfEnergy<Q> operator*(Q alpha) {//multiplication operator overloading
        this->Sigma * alpha;
        return *this;
    }
    SelfEnergy<Q> operator*(double alpha){//multiplication operator overloading
        this->Sigma*alpha;
        return *this;
    }
    SelfEnergy<Q> operator-(const SelfEnergy<Q>& self1){//sum operator overloading
        this->Sigma - self1.Sigma;
        return *this;
    }
};



/*****************************************FUNCTIONS FOR SELF ENERGY*****************************************************/
template <typename Q> Q SelfEnergy<Q>::sval(int iK, int i){
    return Sigma[iK*nSE + i];
}
template <typename Q> Q SelfEnergy<Q>::svalsmooth(int iK, double w){//smoothly interpolates for values between discrete frequency values of mesh

    if(fabs(w)>w_upper_b)
        return 0.;
    else {
        if(fabs(w)!= w_upper_f) {
            int W = fconv_fer(w);
            double x1 = ffreqs[W];
            double x2 = ffreqs[W] + dv;
            double xd = (w - x1) / (x2 - x1);

            Q f1 = sval(iK, W);
            Q f2 = sval(iK, W + 1);

            return (1. - xd) * f1 + xd * f2;
        }
        else if(w == w_upper_f)
            return sval(iK, nSE-1);
        else if(w == w_lower_f)
            return sval(iK, 0);
    }

}
template <typename Q> void SelfEnergy<Q>::setself(int iK, int i, Q val){
    Sigma[iK*nSE + i] = val;
}


#endif //KELDYSH_MFRG_SELFENERGY_H
