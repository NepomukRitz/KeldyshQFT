#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

#include "data_structures.h" // real/complex vector classes
#include <omp.h>             // parallelize initialization of self-energy

/****************** CLASS FOR SELF ENERGY *************/
template <typename Q>
class SelfEnergy{
    // TODO: split into two members: Sigma_R, Sigma_K (?)
    vec<Q> Sigma = vec<Q> (2*nSE*n_in); // factor 2 for Keldysh components: Sigma^R, Sigma^K
public:
    // TODO: comment member functions
    void initialize(Q, Q);
    auto val(int, int, int) const -> Q;
    auto valsmooth(int, double, int) const -> Q;
    void setself(int, int, int, Q);
    void addself(int, int, int, Q);
    auto acc(int) -> Q;// access to the ith element of the vector "SIGMA"
    void direct_set(int, Q);
//operators for self energy

    // TODO: change operator+, operator*
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
template <typename Q> void SelfEnergy<Q>::initialize(Q valR, Q valK) {
    //Assign self energy to initial values
#pragma omp parallel for
    for (int iv=0; iv<nSE; ++iv) {
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setself(0, iv, i_in, valR);
            this->setself(1, iv, i_in, valK);
        }
    }
}

template <typename Q> auto SelfEnergy<Q>::val(int iK, int i, int i_in) const -> Q{
    return Sigma[iK*nSE*n_in + i*n_in + i_in];
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

template <typename Q> auto SelfEnergy<Q>::valsmooth(int iK, double w, int i_in) const -> Q {//smoothly interpolates for values between discrete frequency values of mesh

    if(fabs(w)>w_upper_b)
        //Returns U/2 for Retarded and 0. for Keldysh
        return (1.-(double)iK)*glb_U/2.;
    else {
        if(fabs(w)!= w_upper_f) {
            int W = fconv_fer(w);
            double x1 = ffreqs[W];
            double x2 = ffreqs[W + 1];
            double xd = (w - x1) / (x2 - x1);

            Q f1 = val(iK, W, i_in);
            Q f2 = val(iK, W + 1, i_in);

            return (1. - xd) * f1 + xd * f2;
        }
        else if(w == w_upper_f)
            return val(iK, nSE-1, i_in);
        else if(w == w_lower_f)
            return val(iK, 0, i_in);
    }

}

template <typename Q> void SelfEnergy<Q>::setself(int iK, int i, int i_in, Q val){
    Sigma[iK*nSE + i*n_in + i_in] = val;
}

template <typename Q> void SelfEnergy<Q>::addself(int iK, int i, int i_in, Q val){
    Sigma[iK*nSE + i*n_in + i_in] += val;
}


#endif //KELDYSH_MFRG_SELFENERGY_H
