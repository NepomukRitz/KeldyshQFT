#ifndef KELDYSH_MFRG_SELFENERGY_HPP
#define KELDYSH_MFRG_SELFENERGY_HPP

#include "../../data_structures.hpp" // real/complex vector classes
#include "../../multidimensional/multiarray.hpp"
#include "../n_point/data_buffer.hpp"
#include "../../grids/frequency_grid.hpp"  // interpolate self-energy on new frequency grid
#include "../../utilities/minimizer.hpp"
#include <omp.h>             // parallelize initialization of self-energy
#include "../../symmetries/Keldysh_symmetries.hpp"
#include "../../utilities/write_data2file.hpp"
#include "../../interpolations/InterpolatorLinOrSloppy.hpp"

/// TODO: Use Vertex buffer for self-energy

/****************** CLASS FOR SELF-ENERGY *************/
template <typename Q>
class SelfEnergy{
public:
    using freqGrid_type = bufferFrequencyGrid<selfenergy>;
    using buffer_type = dataBuffer<Q, selfenergy, SE_config.rank, SE_config.num_freqs, SE_config.position_first_freq_index, freqGrid_type, INTERPOLATION>;
    using frequencies_type = std::array<double,1>;
    using index_type = std::array<my_index_t, SE_config.rank>;
 // multidimensional::multiarray<Q,3>;

private:



public:
    buffer_type Sigma;
    Q asymp_val_R = 0.;   //Asymptotic value for the Retarded SE

    explicit  SelfEnergy(double Lambda) : Sigma(Lambda, SE_config.dims) {};
    explicit  SelfEnergy(const freqGrid_type& frequencies_in) : Sigma(0, SE_config.dims) {Sigma.set_VertexFreqGrid( frequencies_in);};

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
    void update_grid(const freqGrid_type&  frequencies_new);   // Interpolate self-energy to updated grid
    /// Interpolate self-energy to input grid with Sigma given by selfEnergy4Sigma
    void update_grid(const freqGrid_type&  frequencies_new, SelfEnergy<Q> selfEnergy4Sigma);
    /// finds optimal grid parameters with minimizer()
    void findBestFreqGrid(bool verbose);       // optimize frequency grid parameters and update self-energy on new grid
    auto shrink_freq_box(double rel_tail_threshold) const -> freqGrid_type ;       // optimize frequency grid parameters and update self-energy on new grid
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
    auto operator+= (Qfac alpha) -> SelfEnergy<Q> {
        this->Sigma += alpha;
        return *this;
    }
    template <typename  Qfac>
    friend SelfEnergy<Q> operator+ (SelfEnergy<Q> lhs, const Qfac& rhs) {
        lhs += rhs;
        return lhs;
    }
    template <typename  Qfac>
    friend SelfEnergy<Q> operator+ (const Qfac& rhs, SelfEnergy<Q> lhs) {
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
    template <typename  Qfac>
    friend SelfEnergy<Q> operator* (const Qfac& rhs, SelfEnergy<Q> lhs) {
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
    /// Elementwise division (needed for error estimate of adaptive ODE solvers)
    auto operator/= (const SelfEnergy<Q>& self1) -> SelfEnergy<Q> {//sum operator overloading
        this->Sigma /= self1.Sigma;
        return *this;
    }
    friend SelfEnergy<Q> operator/ (SelfEnergy<Q> lhs, const SelfEnergy<Q>& rhs) {
        lhs /= rhs;
        return lhs;
    }


    void check_resolution() const;
    void check_symmetries() const;
};



/*****************************************FUNCTIONS FOR SELF-ENERGY*****************************************************/
/**
 * Function initializes the Retarded and Keldysh components of the self-energy
 * @tparam Q    : Type of values (usually comp)
 * @param valR  : Value for the constant Retarded self-energy
 * @param valK  : Value for the constant Keldyh self-energy
 */
template <typename Q> void SelfEnergy<Q>::initialize(Q valR, Q valK) {
    // in particle-hole symmetric_full case (Matsubara formalism) the self-energy vector only stores the imaginary part -> initialize to zero
    // in all other cases: initialize to the Hartree value
    if (KELDYSH || !PARTICLE_HOLE_SYMMETRY){
//#pragma omp parallel for
//        for (int iv=0; iv<nSE; ++iv) {
//            for (int i_in=0; i_in<n_in_K1; ++i_in) {
//                this->setself(0, iv, i_in, valR);
//                if (KELDYSH) this->setself(1, iv, i_in, valK);
//            }
//        }
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
    return Sigma.val(iK, iv, i_in);
}

/**
 * Access function to element i (hdf5-relevant)
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @return Sigma[i]
 */
template <typename Q> auto SelfEnergy<Q>::acc(int i) const -> Q{
    if(i>=0 && i < Sigma.size()) return Sigma.acc(i);
    else {utils::print("Error: Tried to access value outside of self-energy range. Abort."); assert(false);}
}

/**
 * Set value of self-energy at index i
 * @tparam Q : Type of values (usually comp)
 * @param i  : Index
 * @param val: Value
 */
template <typename Q> void SelfEnergy<Q>::direct_set(int i, Q val) {
    if(i>=0 && i < Sigma.get_vec().size()) Sigma.direct_set(i, val);
    else {utils::print("Error: Tried to access value outside of self-energy range. Abort."); assert(false);}
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

    if (std::abs(v) > Sigma.frequencies.get_wupper_b() + inter_tol)    //Check the range of frequency. If too large, return Sigma(\infty)
        //Returns asymptotic value (Hartree contribution for retarded and 0. for Keldysh component)
        return (1.-(double)iK)*(this->asymp_val_R);
    else {
    Q result;
#ifdef DENSEGRID
    result = interpolate_nearest1D<Q>(v, Sigma.frequencies.primary_grid, [&](int i) -> Q {return val(iK, i, i_in);});
#else
    frequencies_type freqs = {v};
        index_type idx;
        idx[my_defs::SE::keldysh]= iK;
        idx[my_defs::SE::internal]= i_in;
    result = Sigma.interpolate_impl(freqs, idx) + (1.-(double)iK)*(this->asymp_val_R);
#endif
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
    Sigma.setvert(val, iK, iv, i_in);
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
     Sigma.setvert(val+Sigma.at(iK, iv, i_in), iK, iv, i_in);
}

template <typename Q> void SelfEnergy<Q>::set_frequency_grid(const SelfEnergy<Q> selfEnergy) {
    Sigma.frequencies = selfEnergy.Sigma.frequencies;
};

template <typename Q> void SelfEnergy<Q>::update_grid(double Lambda) {
    Sigma.update_grid(Lambda);

}

template <typename Q> void SelfEnergy<Q>::update_grid(const freqGrid_type&  frequencies_new) {
    Sigma.update_grid(frequencies_new, Sigma);


}

template <typename Q> void SelfEnergy<Q>::update_grid(const freqGrid_type& frequencies_new, SelfEnergy<Q> selfEnergy4Sigma) {

    Sigma.update_grid(frequencies_new, selfEnergy4Sigma.Sigma);

}



template<typename Q>
class CostSE_Wscale {
    SelfEnergy<Q> selfEnergy_backup;
    bool verbose;
public:
    SelfEnergy<Q> selfEnergy;
    explicit CostSE_Wscale(SelfEnergy<Q> SE_in, bool verbose): selfEnergy(SE_in), selfEnergy_backup(SE_in), verbose(verbose) {
        // remove Hartree contribution
        std::array<size_t, 3> start = {0,0,0};
        std::array<size_t, 3> end   = {0, nFER, n_in};
        auto data = selfEnergy.Sigma.get_vec();
        data.eigen_segment(start, end) -= selfEnergy.asymp_val_R;

        selfEnergy.Sigma.set_vec(data);
        selfEnergy_backup.Sigma.set_vec(data);
        selfEnergy.asymp_val_R = 0.;


        selfEnergy_backup.asymp_val_R = 0.;

    };

    auto operator() (double wscale_test) -> double {
        selfEnergy.Sigma.frequencies.primary_grid.update_Wscale(wscale_test);
        selfEnergy.update_grid(selfEnergy.Sigma.frequencies, selfEnergy_backup);
        double result = selfEnergy.get_curvature_maxSE(verbose);
        /*
        std::string filename = "SE_costCurvature_" + std::to_string(wscale_test) + ".h5";
        rvec v = selfEnergy.frequencies.get_all_frequencies();

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

    Sigma.optimize_grid(verbose);
}

template <typename Q> auto SelfEnergy<Q>::shrink_freq_box(const double rel_tail_threshold) const -> freqGrid_type {
    freqGrid_type frequencies_new = Sigma.shrink_freq_box(rel_tail_threshold, false, 0.);
    return frequencies_new;
}


template <typename Q> double SelfEnergy<Q>::analyze_tails(const bool verbose) const {
    buffer_type Sigma_temp = Sigma;
    //if(KELDYSH) for (int i = 0; i < Sigma.get_vec().size()/2; i++) Sigma_temp.direct_set(i, Sigma_temp.acc(i) - asymp_val_R);
    //else for (int i = 0; i < Sigma.get_vec().size(); i++) Sigma_temp.direct_set(i, Sigma_temp.acc(i) - asymp_val_R);

    double maxabs_SE_total = Sigma_temp.get_vec().max_norm();
    vec<double> maxabsSE_along_w = maxabs(Sigma_temp.get_vec(), Sigma_temp.get_dims(), 1);

    double result = maxabsSE_along_w[0] / maxabs_SE_total;

    if (verbose and mpi_world_rank() == 0)
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
        double max = Sigma.get_vec().max_norm();
        return max;
    }

    else{ //p-norm
        double result = std::abs(Sigma.get_vec().get_elements().pow(p).sum());
        return pow(result, (double)p);
    }
}

/* standard norm: 2-norm */
template <typename Q> auto SelfEnergy<Q>::norm() const -> double {
    return this->norm(0);
}

template <typename Q> auto SelfEnergy<Q>::get_deriv_maxSE(const bool verbose) const -> double {
    //vec<Q> Sigma_temp = Sigma;
    //if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp -= asymp_val_R;
    //else for (int i = 0; i < Sigma.size(); i++) Sigma_temp -= asymp_val_R;

    double max_SE = Sigma.get_deriv_max();
    if (verbose and mpi_world_rank() == 0) {
        std::cout << "max. Derivative in selfenergy:" << std::endl;
        std::cout << "\t  \t" << max_SE << std::endl;
    }

    return max_SE;
}
template <typename Q> auto SelfEnergy<Q>::get_curvature_maxSE(const bool verbose) const -> double {
    //vec<Q> Sigma_temp = Sigma;
    //if(KELDYSH) for (int i = 0; i < Sigma.size()/2; i++) Sigma_temp -= asymp_val_R;
    //else for (int i = 0; i < Sigma.size(); i++) Sigma_temp -= asymp_val_R;

    double max_SE = Sigma.get_curvature_max();
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



template <typename Q> void SelfEnergy<Q>:: check_symmetries() const
{
    if (PARTICLE_HOLE_SYMMETRY) {
        vec<Q> dev(nFER);
        for (int i = 0; i < nFER; i++) {
            for (int j = 0; j < n_in; j++) {
                double v;
                Sigma.frequencies.get_freqs_w(v, i);
                dev[i] = valsmooth(0, v, j) + myconj(valsmooth(0,-v, j)) - 2.*asymp_val_R;
            }
        }

        utils::print("Maximal deviation from PHS: ", dev.max_norm(), "\n");

    }

}

#endif //KELDYSH_MFRG_SELFENERGY_HPP
