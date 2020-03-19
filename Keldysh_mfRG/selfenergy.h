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
    void initialize(Q valR, Q valK);
    auto val(int iK, int iv, int i_in) const -> Q;
    auto valsmooth(int iK, double v, int i_in) const -> Q;
    void setself(int iK, int iv, int i_in, Q val);
    void addself(int iK, int iv, int i_in, Q val);
    auto acc(int i) -> Q;// access to the ith element of the vector "Sigma"
    void direct_set(int i, Q val);

    // operators for self energy
    auto operator+= (const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma += self1.Sigma;
        return *this;
    }
    friend SelfEnergy<Q> operator+ (SelfEnergy<Q> lhs, const SelfEnergy<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (Q alpha) -> SelfEnergy<Q> {
        this->Sigma *= alpha;
        return *this;
    }
    friend SelfEnergy<Q> operator* (SelfEnergy<Q> lhs, const Q& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (double alpha) -> SelfEnergy<Q> {
        this->Sigma *= alpha;
        return *this;
    }
    friend SelfEnergy<Q> operator* (SelfEnergy<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma -= self1.Sigma;
        return *this;
    }
    friend SelfEnergy<Q> operator- (SelfEnergy<Q> lhs, const SelfEnergy<Q>& rhs) {
        lhs -= rhs;
        return lhs;
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

template <typename Q> auto SelfEnergy<Q>::val(int iK, int iv, int i_in) const -> Q{
    return Sigma[iK*nSE*n_in + iv*n_in + i_in];
}

template <typename Q> auto SelfEnergy<Q>::acc(int i) -> Q{
    if(i>=0 && i < Sigma.size()){
    return Sigma[i];}
    else{cout << "Error: Tried to access value outside of self-energy range." << endl;};
}

template <typename Q> void SelfEnergy<Q>::direct_set(int i, Q val) {
    if(i>=0 && i < Sigma.size()){
    Sigma[i] = val;}
    else{cout << "Error: Tried to access value outside of self-energy range." << endl;};
}

template <typename Q> auto SelfEnergy<Q>::valsmooth(int iK, double v, int i_in) const -> Q {//smoothly interpolates for values between discrete frequency values of mesh

    if(abs(v)>w_upper_f)
        //Returns U/2 for retarded and 0. for Keldysh component
        return (1.-(double)iK)*glb_U/2.;
    else {
        if(fabs(v)!= w_upper_f) { // linear interpolation
            int iv = fconv_fer(v); // index corresponding to v
            double x1 = ffreqs[iv]; // lower adjacent frequency value
            double x2 = ffreqs[iv + 1]; // upper adjacent frequency value
            double xd = (v - x1) / (x2 - x1); // distance between adjacent frequnecy values

            Q f1 = val(iK, iv, i_in); // lower adjacent value
            Q f2 = val(iK, iv + 1, i_in);  // upper adjacent value

            return (1. - xd) * f1 + xd * f2; // interpolated value
        }
        else if(v == w_upper_f)
            return val(iK, nSE-1, i_in);
        else if(v == w_lower_f)
            return val(iK, 0, i_in);
    }

}

template <typename Q> void SelfEnergy<Q>::setself(int iK, int iv, int i_in, Q val){
    Sigma[iK*nSE + iv*n_in + i_in] = val;
}

template <typename Q> void SelfEnergy<Q>::addself(int iK, int iv, int i_in, Q val){
    Sigma[iK*nSE + iv*n_in + i_in] += val;
}


#endif //KELDYSH_MFRG_SELFENERGY_H
