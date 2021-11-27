#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

#include "data_structures.h" // real/complex vector classes
#include "grids/frequency_grid.h"  // interpolate self-energy on new frequency grid
#include "minimizer.h"
#include <omp.h>             // parallelize initialization of self-energy
#include "symmetries/Keldysh_symmetries.h"
#include "utilities/write_data2file.h"

/****************** CLASS FOR SELF-ENERGY *************/
template <typename Q>
class SelfEnergy{
    vec<Q> empty_Sigma() {
        return vec<Q> (collapse_all(dimsSE, [](const size_t& a, const size_t& b) -> size_t {return a*b;}));           // only one component in Matsubara formalism
    }

    std::array<size_t,3> dims;

public:
    FrequencyGrid frequencies;
    vec<Q> Sigma = empty_Sigma();
    Q asymp_val_R = 0.;   //Asymptotic value for the Retarded SE

    explicit  SelfEnergy(double Lambda) : frequencies('f', 1, Lambda) {
        for (int i = 0; i < 3; i++) dims[i] = dimsSE[i];};
    explicit  SelfEnergy(const FrequencyGrid& frequencies_in) : frequencies(frequencies_in) {
        for (int i = 0; i < 3; i++) dims[i] = dimsSE[i];};

    void initialize(Q valR, Q valK);    //Initializes SE to given values
    auto val(int iK, int iv, int i_in) const -> Q;  //Returns value at given input on freq grid
    auto valsmooth(int iK, double v, int i_in) const -> Q;  //Returns interpolated value at given input
    void setself(int iK, int iv, int i_in, Q val);  //Sets value of SE at given input
    void addself(int iK, int iv, int i_in, Q val);  //Adds given input to specified location
    auto acc(int i) const -> Q;// access to the ith element of the vector "Sigma" (hdf5-relevant)
    void direct_set(int i, Q val);  //Direct set value to i-th location (hdf5-relevant)
    void set_frequency_grid(const SelfEnergy<Q> selfEnergy);
    auto norm(int p) const -> double;
    auto norm() const -> double;

    /// Interpolate self-energy to updated grid whose grid parameters are multiples of Delta = (Lambda + glb_Gamma)/2
    void update_grid(double Lambda);
    /// Interpolate self-energy to input grid
    void update_grid(FrequencyGrid frequencies_new);   // Interpolate self-energy to updated grid
    /// Interpolate self-energy to input grid with Sigma given by selfEnergy4Sigma
    void update_grid(FrequencyGrid frequencies_new, SelfEnergy<Q> selfEnergy4Sigma);
    /// finds optimal grid parameters with minimizer()
    void findBestFreqGrid(bool verbose);       // optimize frequency grid parameters and update self-energy on new grid
    auto shrink_freq_box(double rel_tail_threshold) const -> FrequencyGrid;       // optimize frequency grid parameters and update self-energy on new grid
    double analyze_tails(bool verbose) const;
    /// computes finite differences of Sigma
    double get_deriv_maxSE(bool verbose) const;
    double get_curvature_maxSE(bool verbose) const;

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

    //FrequencyGrid shrink_freq_box(double rel_tail_threshold) const;
//
    //double analyze_tails(bool verbose) const;

    void check_resolution() const;
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
    if (KELDYSH || !PARTICLE_HOLE_SYMMETRY){
#pragma omp parallel for
        for (int iv=-FREQ_PADDING; iv<nSE+FREQ_PADDING; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                this->setself(0, iv, i_in, valR);
                if (KELDYSH) this->setself(1, iv, i_in, valK);
            }
        }
        this-> asymp_val_R = valR;
    }
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
    if (KELDYSH) return Sigma[iK*(nSE+2*FREQ_PADDING)*n_in + (iv+FREQ_PADDING)*n_in + i_in];
    else         return Sigma[(iv+FREQ_PADDING)*n_in + i_in];
}

/**
 * Access function to element i (hdf5-relevant)
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @return Sigma[i]
 */
template <typename Q> auto SelfEnergy<Q>::acc(int i) const -> Q{
    if(i>=0 && i < Sigma.size()) return Sigma[i];
    else {print("Error: Tried to access value outside of self-energy range. Abort."); assert(false);}
}

/**
 * Set value of self-energy at index i
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @param val: Value
 */
template <typename Q> void SelfEnergy<Q>::direct_set(int i, Q val) {
    if(i>=0 && i < Sigma.size()) Sigma[i] = val;
    else {print("Error: Tried to access value outside of self-energy range. Abort."); assert(false);}
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
    Q result;
    if (INTERPOLATION == linear) result = interpolate_lin1D<Q>(v, frequencies, [&](int i) -> Q {return val(iK, i, i_in);});
    else result = interpolate_lin_on_aux1D<Q>(v, frequencies, [&](int i) -> Q {return val(iK, i, i_in);});
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
    Sigma[iK*(nSE+2*FREQ_PADDING) + (iv+FREQ_PADDING)*n_in + i_in] = val;
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
    Sigma[iK*(nSE+2*FREQ_PADDING) + (iv+FREQ_PADDING)*n_in + i_in] += val;
}

template <typename Q> void SelfEnergy<Q>::set_frequency_grid(const SelfEnergy<Q> selfEnergy) {
    this->frequencies = selfEnergy.frequencies;
};

template <typename Q> void SelfEnergy<Q>::update_grid(double Lambda) {
    FrequencyGrid frequencies_new = this->frequencies; // new frequency grid
    frequencies_new.rescale_grid(Lambda);              // rescale new frequency grid

    vec<Q> Sigma_new = empty_Sigma();                     // temporary self-energy vector
    for (int iK=0; iK<nK_SE; ++iK) {
        if (!KELDYSH && (iK == 1)) break; // Only Keldysh index 0 for Matsubara
        for (int iv=0; iv<nSE; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                // interpolate old values to new vector
                Sigma_new[iK*(nSE+2*FREQ_PADDING)*n_in + (iv*FREQ_PADDING)*n_in + i_in] = this->valsmooth(iK, frequencies_new.get_ws(iv), i_in);
            }
        }
#if FREQ_PADDING == 1
        for (int i_in=0; i_in<n_in; ++i_in) {
            // set asymptotic values
            Sigma_new[0*nSE*n_in + 0*n_in + i_in] = asymp_val_R;
            Sigma_new[0*nSE*n_in + (nSE+FREQ_PADDING)*n_in + i_in] = asymp_val_R;
        }
#endif
    }
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}

template <typename Q> void SelfEnergy<Q>::update_grid(FrequencyGrid frequencies_new) {

    vec<Q> Sigma_new (nK_SE*(nSE+2*FREQ_PADDING)*n_in);                     // temporary self-energy vector
    for (int iK=0; iK<nK_SE; ++iK) {
        if (!KELDYSH && (iK == 1)) break; // Only Keldysh index 0 for Matsubara
        for (int iv=0; iv<nSE; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                // interpolate old values to new vector
                Sigma_new[iK*(nSE+2*FREQ_PADDING)*n_in + (iv+FREQ_PADDING)*n_in + i_in] = this->valsmooth(iK, frequencies_new.get_ws(iv), i_in);
            }
        }
#if FREQ_PADDING == 1
        for (int i_in=0; i_in<n_in; ++i_in) {
            // set asymptotic values
            Sigma_new[0*nSE*n_in + 0*n_in + i_in] = asymp_val_R;
            Sigma_new[0*nSE*n_in + (nSE+FREQ_PADDING)*n_in + i_in] = asymp_val_R;
        }
#endif
    }
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}

template <typename Q> void SelfEnergy<Q>::update_grid(FrequencyGrid frequencies_new, SelfEnergy<Q> selfEnergy4Sigma) {

    vec<Q> Sigma_new (nK_SE*(nSE+2*FREQ_PADDING)*n_in);                     // temporary self-energy vector

    for (int iK=0; iK<nK_SE; ++iK) {
        if (!KELDYSH && (iK == 1)) break; // Only Keldysh index 0 for Matsubara
        for (int iv=0; iv<nSE; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                // interpolate old values to new vector
                Sigma_new[iK*(nSE+2*FREQ_PADDING)*n_in + (iv+FREQ_PADDING)*n_in + i_in] = selfEnergy4Sigma.valsmooth(iK, frequencies_new.get_ws(iv), i_in);
            }
        }
#if FREQ_PADDING == 1
        for (int i_in=0; i_in<n_in; ++i_in) {
            // set asymptotic values
            Sigma_new[0*nSE*n_in + 0*n_in + i_in] = asymp_val_R;
            Sigma_new[0*nSE*n_in + (nSE+FREQ_PADDING)*n_in + i_in] = asymp_val_R;
        }
#endif
    }
    this->frequencies = frequencies_new; // update frequency grid to new rescaled grid
    this->Sigma = Sigma_new;             // update selfenergy to new interpolated values
}



template<typename Q>
class CostSE_wupper {
    SelfEnergy<Q> selfEnergy;
    double rel_tailsize = 1e-2;
public:
    explicit CostSE_wupper(SelfEnergy<Q> SE_in): selfEnergy(SE_in) {
        // remove Hartree contribution
        for (int iv=0; iv<nSE+2*FREQ_PADDING; ++iv) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                selfEnergy.Sigma[iv*n_in + i_in] -= selfEnergy.asymp_val_R;
            }
        }
        selfEnergy.asymp_val_R = 0.;
    };

    auto operator() (double w_upper_test) -> double {

        double max = selfEnergy.norm(0);
        if (w_upper_test < selfEnergy.frequencies.w_upper) {
            double result = std::abs((std::abs(selfEnergy.valsmooth(0, w_upper_test, 0)) +
                                      std::abs(selfEnergy.valsmooth(0, -w_upper_test, 0))) / max - rel_tailsize);
            return result;
        }
        else {
            double tupper_test = selfEnergy.frequencies.grid_transf(w_upper_test);
            double result_tmp = std::abs(selfEnergy.Sigma[0]) + std::abs(selfEnergy.Sigma[nFER-1]);
            double factor = ((1. - tupper_test)/(1.-selfEnergy.frequencies.t_upper));
            double result = std::abs(result_tmp * factor / max - rel_tailsize);
            return result;
        }
    }
};

template<typename Q>
class CostSE_Wscale {
    SelfEnergy<Q> selfEnergy_backup;
    bool verbose;
public:
    SelfEnergy<Q> selfEnergy;
    explicit CostSE_Wscale(SelfEnergy<Q> SE_in, bool verbose): selfEnergy(SE_in), selfEnergy_backup(SE_in), verbose(verbose) {
        // remove Hartree contribution
        for (int iv = 0; iv < nSE + 2*FREQ_PADDING; ++iv) {
            for (int i_in = 0; i_in < n_in; ++i_in) {
                selfEnergy.Sigma[iv * n_in + i_in] -= selfEnergy.asymp_val_R;
            }
        }

        selfEnergy.asymp_val_R = 0.;

        for (int iv = 0; iv < nSE + 2*FREQ_PADDING; ++iv) {
            for (int i_in = 0; i_in < n_in; ++i_in) {
                selfEnergy_backup.Sigma[iv * n_in + i_in] -= selfEnergy_backup.asymp_val_R;
            }
        }

        selfEnergy_backup.asymp_val_R = 0.;

    };

    auto operator() (double wscale_test) -> double {
        selfEnergy.frequencies.update_Wscale(wscale_test);
        selfEnergy.update_grid(selfEnergy.frequencies, selfEnergy_backup);
        double result = selfEnergy.get_curvature_maxSE(verbose);
        /*
        std::string filename = "SE_costCurvature_" + std::to_string(wscale_test) + ".h5";
        rvec v = selfEnergy.frequencies.get_ws_vec();

        rvec SE_re = selfEnergy.Sigma.real();
        rvec SE_im = selfEnergy.Sigma.imag();
        write_h5_rvecs(filename,
                       {"v", "SE_re", "SE_im"},
                       {v, SE_re, SE_im});
        */

        return result;



    }

};

template <typename Q> void SelfEnergy<Q>::findBestFreqGrid(const bool verbose) {
    double rel_tail_threshold = 1e-2;

    FrequencyGrid frequencies_new = shrink_freq_box(rel_tail_threshold);
    update_grid(frequencies_new);

    SelfEnergy<Q> SEtemp = *this;
    //SEtemp.update_grid(Lambda);

    //double wmax_current = SEtemp.frequencies.w_upper;
    double a_Wscale = SEtemp.frequencies.W_scale / 10.;
    double m_Wscale = SEtemp.frequencies.W_scale;
    double b_Wscale = SEtemp.frequencies.W_scale * 10;
    CostSE_Wscale<Q> cost(SEtemp, verbose);
    minimizer(cost, a_Wscale, m_Wscale, b_Wscale, 100, verbose, false, 1., 0.);
    frequencies_new.update_Wscale(m_Wscale);

    update_grid(frequencies_new);


}

template <typename Q> auto SelfEnergy<Q>::shrink_freq_box(const double rel_tail_threshold) const -> FrequencyGrid {
    vec<Q> Sigma_temp = Sigma;
    if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp[i] -= asymp_val_R;
    else for (int i = 0; i < Sigma.size(); i++) Sigma_temp[i] -= asymp_val_R;

    double maxmax = Sigma.max_norm();
    vec<double> maxabsSE_along_w = maxabs(Sigma_temp, dims, 1) * (1/maxmax);

    FrequencyGrid frequencies_new = freqGrid::shrink_freq_box(frequencies, rel_tail_threshold, maxabsSE_along_w);

    return frequencies_new;
}


template <typename Q> double SelfEnergy<Q>::analyze_tails(const bool verbose) const {
    vec<Q> Sigma_temp = Sigma;
    if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp[i] -= asymp_val_R;
    else for (int i = 0; i < Sigma.size(); i++) Sigma_temp[i] -= asymp_val_R;

    double maxabs_SE_total = Sigma_temp.max_norm();
    vec<double> maxabsSE_along_w = maxabs(Sigma_temp, dims, 1);

    double result = maxabsSE_along_w[FREQ_PADDING] / maxabs_SE_total;

    if (verbose)
    {
        std::cout << "rel. magnitude of tails in self energy :";
        std::cout << "\t\t" << result << std::endl;
    }

    return result;
}
/*
 * p-norm for the SelfEnergy
 */
template <typename Q> auto SelfEnergy<Q>::norm(const int p) const -> double {
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
template <typename Q> auto SelfEnergy<Q>::norm() const -> double {
    return this->norm(2);
}

template <typename Q> auto SelfEnergy<Q>::get_deriv_maxSE(const bool verbose) const -> double {
    //vec<Q> Sigma_temp = Sigma;
    //if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp -= asymp_val_R;
    //else for (int i = 0; i < Sigma.size(); i++) Sigma_temp -= asymp_val_R;

    double maxmax = Sigma.max_norm();
    double dt = frequencies.dt;

    double max_SE = (::power2(::partial_deriv<Q,3>(Sigma, frequencies.get_ts_vec(), dims, 1)*dt*(1/maxmax))).max_norm();

    if (verbose) {
        std::cout << "max. Derivative in selfenergy:" << std::endl;
        std::cout << "\t  \t" << max_SE << std::endl;
    }

    return max_SE;
}
template <typename Q> auto SelfEnergy<Q>::get_curvature_maxSE(const bool verbose) const -> double {
    //vec<Q> Sigma_temp = Sigma;
    //if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp -= asymp_val_R;
    //else for (int i = 0; i < Sigma.size(); i++) Sigma_temp -= asymp_val_R;

    double maxmax = Sigma.max_norm();
    double dt = frequencies.dt;
    //double max_SE = ::power2(::get_finite_differences(Sigma)).max_norm();
    //return max_SE;
    const std::array<size_t,3> dims1 = {n_in, nK_SE, nFER+2*FREQ_PADDING};
    const std::array<size_t,3> perm1 = {2, 0, 1};
    double max_SE = (::power2(::partial_deriv<Q,3>(::partial_deriv<Q,3>(Sigma, frequencies.get_ts_vec(), dims, 1), frequencies.get_ts_vec(), dims, 1)*dt*dt*(1/maxmax))).max_norm();

    if (verbose and mpi_world_rank() == 0) {
        std::cout << "max. Curvature in SE:";
        std::cout << "\t  \t" << max_SE << std::endl;
    }

    return max_SE;


}


template <typename Q> void SelfEnergy<Q>:: check_resolution() const
{
    double derivmax_SE = get_deriv_maxSE(true);
    double curvmax_SE = get_curvature_maxSE(true);

    /// TODO: dump state in file if certain thresholds are exceeded
}

#endif //KELDYSH_MFRG_SELFENERGY_H
