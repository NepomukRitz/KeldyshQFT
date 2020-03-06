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
    auto sval(int, int) -> Q;
    auto svalsmooth(int, double) -> Q;
    void setself(int, int, Q);
    void addself(int, int, Q);
    auto acc(int) ->Q;// access to the ith element of the vector "SIGMA"
    void direct_set(int,Q);
//operators for self energy

    auto operator+(const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma + self1.Sigma;
        return *this;
    }
    auto operator+=(const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma += self1.Sigma;
        return *this;
    }
    auto operator*(Q alpha) -> SelfEnergy<Q> {//multiplication operator overloading
        this->Sigma * alpha;
        return *this;
    }
    auto operator*(double alpha) -> SelfEnergy<Q> {
        this->Sigma*alpha;
        return *this;
    }
    auto operator*=(double alpha) -> SelfEnergy<Q> {
        this->Sigma*=alpha;
        return *this;
    }
    auto operator-(const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma - self1.Sigma;
        return *this;
    }
    auto operator-=(const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma -= self1.Sigma;
        return *this;
    }

};



/*****************************************FUNCTIONS FOR SELF ENERGY*****************************************************/
template <typename Q> auto SelfEnergy<Q>::sval(int iK, int i) -> Q{
    return Sigma[iK*nSE + i];
}

template <typename Q> auto SelfEnergy<Q>::acc(int i) -> Q{
    if(i>=0 && i < Sigma.size()){
    return Sigma[i];}
    else{cout << "Error: Tried to access value outside of self energy range" << endl;};
}

template <typename Q> void SelfEnergy<Q>::direct_set(int i, Q value) {
    if(i>=0 && i < Sigma.size()){
    Sigma[i] = value;}
    else{cout << "Error: Tried to access value outside of self energy range" << endl;};
}
template <typename Q> auto SelfEnergy<Q>::svalsmooth(int iK, double w) -> Q{//smoothly interpolates for values between discrete frequency values of mesh

    if(fabs(w)>w_upper_b)
        //Returns U/2 for Retarded and 0. for Keldysh
        return (1.-(double)iK)*glb_U/2.;
    else {
        if(fabs(w)!= w_upper_f) {
            int W = fconv_fer(w);
            double x1 = ffreqs[W];
            double x2 = ffreqs[W + 1];
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
template <typename Q> void SelfEnergy<Q>::addself(int iK, int i, Q val){
    Sigma[iK*nSE + i] += val;
}


#endif //KELDYSH_MFRG_SELFENERGY_H
