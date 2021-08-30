#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

#include "data_structures.h" // real/complex vector classes
#include "grids/frequency_grid.h"  // interpolate self-energy on new frequency grid
#include "minimizer.h"
#include <omp.h>             // parallelize initialization of self-energy

/****************** CLASS FOR SELF-ENERGY *************/
template <typename Q>
class SelfEnergy{
public:
    FrequencyGrid frequencies;
#ifdef KELDYSH_FORMALISM
    vec<Q> Sigma = vec<Q> (2*nSE*n_in); // factor 2 for Keldysh components: Sigma^R, Sigma^K
#else
    vec<Q> Sigma = vec<Q> (nSE*n_in); // only one component in Matsubara formalism
#endif
    Q asymp_val_R = 0.;   //Asymptotic value for the Retarded SE

    explicit  SelfEnergy(double Lambda) : frequencies('f', 1, Lambda) {};
    explicit  SelfEnergy(const FrequencyGrid& frequencies_in) : frequencies(frequencies_in) {};


    void initialize(Q valR, Q valK);    //Initializes SE to given values
    auto val(int iK, int iv, int i_in) const -> Q;  //Returns value at given input on freq grid
    auto valsmooth(int iK, double v, int i_in) const -> Q;  //Returns interpolated value at given input
    void setself(int iK, int iv, int i_in, Q val);  //Sets value of SE at given input
    void addself(int iK, int iv, int i_in, Q val);  //Adds given input to specified location
    auto acc(int i) -> Q;// access to the ith element of the vector "Sigma" (hdf5-relevant)
    void direct_set(int i, Q val);  //Direct set value to i-th location (hdf5-relevant)
    void set_frequency_grid(const SelfEnergy<Q>& selfEnergy);
    void update_grid(double Lambda);  // Interpolate self-energy to updated grid
    void update_grid(FrequencyGrid frequencies_new);   // Interpolate self-energy to updated grid
    void findBestFreqGrid(double Lambda);       // optimize frequency grid parameters and update self-energy on new grid
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

    double get_deriv_maxSE() const;

    double cost_wupper(double w_upper_test, void *params);
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
    else{std::cout << "Error: Tried to access value outside of self-energy range." << std::endl;};
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
    else{std::cout << "Error: Tried to access value outside of self-energy range." << std::endl;};
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

    if (std::abs(v) > this->frequencies.w_upper)    //Check the range of frequency. If too large, return Sigma(\infty)
        //Returns asymptotic value (Hartree contribution for retarded and 0. for Keldysh component)
        return (1.-(double)iK)*(this->asymp_val_R);
    else {
            Q result = interpolate1D<Q>(v, frequencies, [&](int i) -> Q {return val(iK, i, i_in);});
            return result;
    }

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
                Sigma_new[iK*nSE*n_in + iv*n_in + i_in] = this->valsmooth(iK, frequencies_new.ws[iv], i_in);
            }
        }
#ifdef KELDYSH_FORMALISM
    }
#endif
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}

template <typename Q> void SelfEnergy<Q>::update_grid(FrequencyGrid frequencies_new) {

    vec<Q> Sigma_new (2*nSE*n_in);                     // temporary self-energy vector
#ifdef KELDYSH_FORMALISM
    for (int iK=0; iK<2; ++iK) {
#else
        int iK = 0;
#endif
        for (int iv=0; iv<nSE; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                // interpolate old values to new vector
                Sigma_new[iK*nSE*n_in + iv*n_in + i_in] = this->valsmooth(iK, frequencies_new.ws[iv], i_in);
            }
        }
#ifdef KELDYSH_FORMALISM
    }
#endif
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}



template<typename Q>
class Cost_wupper {
    SelfEnergy<Q> selfEnergy;
    double rel_tailsize = 1e-3;
public:
    explicit Cost_wupper(SelfEnergy<Q> SE_in): selfEnergy(SE_in) {};

    auto operator() (double w_upper_test) -> double {
        {
            double max = selfEnergy.norm(0);
            if (w_upper_test < selfEnergy.frequencies.w_upper) {
                double result = std::abs((std::abs(selfEnergy.valsmooth(0, w_upper_test, 0)) +
                                          std::abs(selfEnergy.valsmooth(0, -w_upper_test, 0))) / max - rel_tailsize);
                return result;
            }
            else {
                double tupper_test = selfEnergy.frequencies.grid_transf(w_upper_test);
                return std::abs(std::abs(selfEnergy.Sigma[0]) + std::abs(selfEnergy.Sigma[nFER-1]) * (1. - tupper_test/selfEnergy.frequencies.t_upper) / max - rel_tailsize);
            }


        }
    }
};
template <typename Q> void SelfEnergy<Q>::findBestFreqGrid(double Lambda) {

    SelfEnergy<Q> SEtemp = *this;
    SEtemp.update_grid(Lambda);


    double a_wupper = 0.;
    double m_wupper = SEtemp.frequencies.w_upper;
    double b_wupper = SEtemp.frequencies.w_upper * 100;
    Cost_wupper<Q> cost(SEtemp);
    minimizer(cost, a_wupper, m_wupper, b_wupper, 100);

    update_grid(SEtemp.frequencies);


}


/*
 * p-norm for the SelfEnergy
 */
template <typename Q> auto SelfEnergy<Q>::norm(const int p) -> double {
    if(p==0){ //max norm
        double max = 0.;
        for (auto value : (this->Sigma)){
            if(std::abs(value) > max){
                max = std::abs(value);
            }
        }
        return max;
    }

    else{ //p-norm
        double result = 0;
        for (auto value : (this->Sigma)){
            result += pow(std::abs(value), (double)p);
        }
        return pow(result, 1./((double)p));
    }
}

/* standard norm: 2-norm */
template <typename Q> auto SelfEnergy<Q>::norm() -> double {
    return this->norm(2);
}

template <typename Q> auto SelfEnergy<Q>::get_deriv_maxSE() const -> double {
    double max_SE = ::power2(::get_finite_differences(Sigma)).max_norm();
    return max_SE;

}

#endif //KELDYSH_MFRG_SELFENERGY_H
