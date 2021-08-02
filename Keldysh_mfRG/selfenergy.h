#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

#include "data_structures.h" // real/complex vector classes
#include "frequency_grid.h"  // interpolate self-energy on new frequency grid
#include <omp.h>             // parallelize initialization of self-energy

/****************** CLASS FOR SELF-ENERGY *************/
template <typename Q>
class SelfEnergy{
public:
    FrequencyGrid frequencies;
#ifdef KELDYSH_FORMALISM
    // TODO: split into two members: Sigma_R, Sigma_K (?)
    vec<Q> Sigma = vec<Q> (2*nSE*n_in); // factor 2 for Keldysh components: Sigma^R, Sigma^K
#else
    vec<Q> Sigma = vec<Q> (nSE*n_in); // only one component in Matsubara formalism
#endif
    Q asymp_val_R = 0.;   //Asymptotic value for the Retarded SE

    SelfEnergy() : frequencies('f', 1) {};
    SelfEnergy(double Lambda) : frequencies('f', 1, Lambda) {};

    void initialize(Q valR, Q valK);    //Initializes SE to given values
    auto val(int iK, int iv, int i_in) const -> Q;  //Returns value at given input on freq grid
    auto valsmooth(int iK, double v, int i_in) const -> Q;  //Returns interpolated value at given input
    void setself(int iK, int iv, int i_in, Q val);  //Sets value of SE at given input
    void addself(int iK, int iv, int i_in, Q val);  //Adds given input to specified location
    auto acc(int i) -> Q;// access to the ith element of the vector "Sigma" (hdf5-relevant)
    void direct_set(int i, Q val);  //Direct set value to i-th location (hdf5-relevant)
    void set_frequency_grid(const SelfEnergy<Q>& selfEnergy);
    void update_grid(double Lambda);  // Interpolate self-energy to updated grid
    auto norm(int p) -> double;
    auto norm() -> double;

    // operators for self-energy
    auto operator+= (const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma += self1.Sigma;
        return *this;
    }
    friend SelfEnergy<Q> operator+ (SelfEnergy<Q> lhs, const SelfEnergy<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    template <typename  Qfac>
    auto operator*= (Qfac alpha) -> SelfEnergy<Q> {
        this->Sigma *= alpha;
        return *this;
    }
    template <typename  Qfac>
    friend SelfEnergy<Q> operator* (SelfEnergy<Q> lhs, const Qfac& rhs) {
        lhs *= rhs;
        return lhs;
    }
    //auto operator*= (double alpha) -> SelfEnergy<Q> {
    //    this->Sigma *= alpha;
    //    return *this;
    //}
    //friend SelfEnergy<Q> operator* (SelfEnergy<Q> lhs, const double& rhs) {
    //    lhs *= rhs;
    //    return lhs;
   // }
    auto operator-= (const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma -= self1.Sigma;
        return *this;
    }
    friend SelfEnergy<Q> operator- (SelfEnergy<Q> lhs, const SelfEnergy<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }

};



/*****************************************FUNCTIONS FOR SELF-ENERGY*****************************************************/
/**
 * Function initializes the Retarded and Keldysh components of the self-energy
 * @tparam Q    : Type of values (usually comp)
 * @param valR  : Value for the constant Retarded self-energy
 * @param valK  : Value for the constant Keldyh self-energy
 */
template <typename Q> void SelfEnergy<Q>::initialize(Q valR, Q valK) {
// in particle-hole symmetric case (Matsubara formalism) the self-energy vector only stores the imaginary part -> initialize to zero
// in all other cases: initialize to the Hartree value
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
#pragma omp parallel for

    for (int iv=0; iv<nSE; ++iv) {
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setself(0, iv, i_in, valR);
#ifdef KELDYSH_FORMALISM
            this->setself(1, iv, i_in, valK);
#endif
        }
    }
    this-> asymp_val_R = valR;
#endif
}

/**
 * Returns the saved value of the self-energy
 * @tparam Q    : Type of values (usually comp)
 * @param iK    : Keldysh index (either 0 or 1)
 * @param iv    : Frequency index
 * @param i_in  : Internal structure index
 * @return The value of SigmaR/K at the chosen indices
 */
template <typename Q> auto SelfEnergy<Q>::val(int iK, int iv, int i_in) const -> Q{
#ifdef KELDYSH_FORMALISM
    return Sigma[iK*nSE*n_in + iv*n_in + i_in];
#else
    return Sigma[iv*n_in + i_in];
#endif
}

/**
 * Access function to element i (hdf5-relevant)
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @return Sigma[i]
 */
template <typename Q> auto SelfEnergy<Q>::acc(int i) -> Q{
    if(i>=0 && i < Sigma.size()){
    return Sigma[i];}
    else{cout << "Error: Tried to access value outside of self-energy range." << endl;};
}

/**
 * Set value of self-energy at index i
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @param val: Value
 */
template <typename Q> void SelfEnergy<Q>::direct_set(int i, Q val) {
    if(i>=0 && i < Sigma.size()){
    Sigma[i] = val;}
    else{cout << "Error: Tried to access value outside of self-energy range." << endl;};
}

/**
 * Linear interpolating function returning the value of the Retarded or Keldysh self-energy at freq v
 * @tparam Q    : Type of values (usually comp)
 * @param iK    : Keldysh index (either 0 or 1)
 * @param v     : Frequency
 * @param i_in  : Internal index
 * @return      : Interpolated value using the saved values as basis
 */
template <typename Q> auto SelfEnergy<Q>::valsmooth(int iK, double v, int i_in) const -> Q {//smoothly interpolates for values between discrete frequency values of mesh
    double small = 1e-12;
    double tv = this->frequencies.grid_transf(v);
    if (fabs(tv) < this->frequencies.W_upper + small) {    //Check the range of frequency. If too large, return Sigma(\infty)
                                                            //Returns U/2 for retarded and 0. for Keldysh component
        //if (fabs(v) != this->frequencies.w_upper) {
        int iv = (int) ((tv - this->frequencies.W_lower) / this->frequencies.dW); // index corresponding to v

        double x1 = this->frequencies.Ws[iv]; // lower adjacent frequency value
        double x2 = this->frequencies.Ws[iv + 1]; // upper adjacent frequency value
        double xd = (tv - x1) / (x2 - x1); // distance between adjacent frequnecy values

        Q f1 = val(iK, iv, i_in); // lower adjacent value
        Q f2 = val(iK, iv + 1, i_in);  // upper adjacent value

        Q result = (1. - xd) * f1 + xd * f2;
        assert(isfinite(result));
        return result; // interpolated value
        //}
        //else if (v == this->frequencies.w_upper) //Exactly at the upper boundary
        //    return val(iK, nSE-1, i_in);
        //else if (v == this->frequencies.w_lower) //Exactly at the lower boundary
        //    return val(iK, 0, i_in);
    }
    else
        return (1.-(double)iK)*(this->asymp_val_R);
}

/**
 * Sets the value of the self-energy to a specific value at the specified location
 * @tparam Q    : Type of values (usually comp)
 * @param iK    : Keldysh index (either 0 or 1)
 * @param iv    : Frequency index
 * @param i_in  : Internal index
 * @param val   : Value of the self-energy at this given point
 */
template <typename Q> void SelfEnergy<Q>::setself(int iK, int iv, int i_in, Q val){
    Sigma[iK*nSE + iv*n_in + i_in] = val;
}

/**
 * Adds the given value at the specified location
 * @tparam Q    : Type of values (usually comp)
 * @param iK    : Keldysh index (either 0 or 1)
 * @param iv    : Frequency index
 * @param i_in  : Internal index
 * @param val   : Value of the self-energy at this given point
 */
template <typename Q> void SelfEnergy<Q>::addself(int iK, int iv, int i_in, Q val){
    Sigma[iK*nSE + iv*n_in + i_in] += val;
}

template <typename Q> void SelfEnergy<Q>::set_frequency_grid(const SelfEnergy<Q>& selfEnergy) {
    this->frequencies = selfEnergy.frequencies;
};

template <typename Q> void SelfEnergy<Q>::update_grid(double Lambda) {
    FrequencyGrid frequencies_new = this->frequencies; // new frequency grid
    frequencies_new.rescale_grid(Lambda);              // rescale new frequency grid

    vec<Q> Sigma_new (2*nSE*n_in);                     // temporary self-energy vector
#ifdef KELDYSH_FORMALISM
    for (int iK=0; iK<2; ++iK) {
#else
    int iK = 0;
#endif
        for (int iv=0; iv<nSE; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                // interpolate old values to new vector
                Sigma_new[iK*nSE*n_in + iv*n_in + i_in] = this->valsmooth(iK, frequencies_new.w[iv], i_in);
            }
        }
#ifdef KELDYSH_FORMALISM
    }
#endif
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}

/*
 * p-norm for the SelfEnergy
 */
template <typename Q> auto SelfEnergy<Q>::norm(const int p) -> double {
    if(p==0){ //max norm
        double max = 0.;
        for (auto value : (this->Sigma)){
            if(abs(value) > max){
                max = abs(value);
            }
        }
        return max;
    }

    else{ //p-norm
        double result = 0;
        for (auto value : (this->Sigma)){
            result += pow(abs(value), (double)p);
        }
        return pow(result, 1./((double)p));
    }
}

/* standard norm: 2-norm */
template <typename Q> auto SelfEnergy<Q>::norm() -> double {
    return this->norm(2);
}


#endif //KELDYSH_MFRG_SELFENERGY_H
