/**
 * Set up the frequency grid //
 * Functions that initialize the grid and provide conversion between doubles and grid indices.
 * The grid points are obtained as follows:
 *      The values of t are evenly distributed between t_lower and t_upper.
 *      The frequencies w(t) are obtained according to a function.
 * Three different grid types:
 * GRID=0: non-linear grid according to a function in every direction:
 *  W_scale determines 'how non-linear' the grid behaves.
 *     version 1: w(t) = W_scale t  /sqrt(1 - t^2)
 *     version 2: w(t) = W_scale t^2/sqrt(1 - t^2)
 *     version 3: w(t) = W_scale t  /    (1 - t^2)
 *     version 4: w(t) = W_scale t^2/    (1 - t^2)
 * GRID=1: hybrid grid with
 *     quadratic part -- linear part -- rational part
 * GRID=2: polar coordinates
 *     2D: (w/2, v) = rho * (cos phi, sin phi)
 *     3D: /// TODO
 */

/**
 * Things to do when introducing a new grid:
 *      write new template specialization for FrequencyGrid
 *      adapt hdf5_utilities
 *      adapt grid optimization
 */

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_HPP
#define KELDYSH_MFRG_FREQUENCY_GRID_HPP

#include <cmath>        // for sqrt, log, exp
#include "../data_structures.hpp"
#include "../utilities/util.hpp"
#include "../parameters/master_parameters.hpp" // for frequency/Lambda limits and number of frequency/Lambda points
#include "../utilities/math_utils.hpp"
#include "../symmetries/Keldysh_symmetries.hpp"
#include <cassert>
#include "H5Cpp.h"

// TODO(low): implement functions used for GRID=3 also for GRID=1,2,4
/// TODO: implement hybrid frequency grid, improve/unify treatment of different frequency meshes and their frequency parameters

template<K_class k, typename Q> class vertexDataContainer; // forward declaration
template<typename Q, bool differentiated> class State; // forward declaration

#define PARAMETRIZED_GRID
#if not KELDYSH_FORMALISM and not ZERO_TEMP
#define DENSEGRID
#endif
#ifdef DENSEGRID
const bool dense = true;
#else
const bool dense = false;
#endif

#define HYBRID_GRID_OPTION 0


namespace hdf5_impl {
    template<typename gridType> void init_freqgrid_from_hdf_LambdaLayer(H5::Group &group, gridType &freqgrid, const int Lambda_it, const double Lambda);
    template<typename gridType> void write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const gridType& freqgrid, const int Lambda_it, const int numberLambdaLayers, const bool file_exists, const bool verbose);
}

enum frequencyGridType {eliasGrid, hybridGrid, angularGrid};

template <frequencyGridType freqGridType>
class FrequencyGrid {};

template <>
class FrequencyGrid<eliasGrid> {
    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type> friend class DataContainer;
    template<typename gridType> friend void hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, gridType& freqgrid, int Lambda_it, double Lambda);
    template<typename gridType> friend void hdf5_impl::write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const gridType& freqgrid, const int Lambda_it, const int numberLambdaLayers, const bool file_exists, const bool verbose);



public:
    /// essential grid parameters:
    int number_of_gridpoints;
    freqType w_upper;             // lower bound of frequency grid
    freqType w_lower;             // lower bound of frequency grid
    freqType W_scale;             // non-linearity of t_from_frequency()
    freqType w_center = 0;        // symmetry point of frequency grid

    /// guess essential parameters from value of Lambda
    void guess_essential_parameters(double Lambda, const fRG_config& config);


//private:
    /// grid identifier
    char type;
    unsigned int diag_class;
    bool purely_positive;

    /// auxiliary grid parameters:
    freqType t_upper;                     // upper bound of auxiliary grid
    freqType t_lower;                     // lower bound of auxiliary grid
    freqType spacing_auxiliary_gridpoint; // spacing on linear auxiliary grid
    double U_factor = 0./3.;   // determines scale_factor()
    double Delta_factor = 10.;  // determines scale_factor()
    void derive_auxiliary_parameters();

    /// list of all frequencies:
    vec<freqType> all_frequencies;                     // frequency grid
    vec<freqType> auxiliary_grid;                    // linear auxiliary grid (related to all_frequencies by all_frequencies(auxiliary_grid)=frequency_from_t(auxiliary_grid))



public:

    ///This constructor initializes a frequency grid with the global values. This is not needed anymore!
    FrequencyGrid(char type_in, unsigned int diag_class_in, double Lambda, const fRG_config& config, bool purely_positive=false) : type(type_in), diag_class(diag_class_in), purely_positive(purely_positive) {
        guess_essential_parameters(Lambda, config);
    };

    /// getter functions
    auto get_frequency(int index) const -> freqType {assert(index>=0); assert(index<number_of_gridpoints); assert(isfinite(all_frequencies[index])); return all_frequencies[index];};
    auto get_auxiliary_gridpoint(int index) const -> freqType {assert(index>=0); assert(index<number_of_gridpoints); assert(isfinite(auxiliary_grid[index])); return auxiliary_grid[index];};
    auto get_all_frequencies() const -> vec<freqType> {return all_frequencies;}
    auto get_all_auxiliary_gridpoints() const -> vec<freqType> {return auxiliary_grid;}
    freqType get_spacing_auxiliary_gridpoints() const {return spacing_auxiliary_gridpoint;}
    char get_type() const {return type;};
    char get_diag_class() const {return diag_class;};

    /// setter functions
    // set w_upper and W_scale with some checks for Matsubara T>0:
    void set_essential_parameters(const freqType wmax_in, const freqType Wscale_in) {

        if constexpr (!KELDYSH && !ZERO_T){
#ifndef DENSEGRID
            /// TODO: Figure out better solution for non-uniform grid
            // for Matsubara T>0: pick grid such that no frequencies occur twice
            if (type == 'b') {
                w_upper = std::max(round2bfreq(wmax_in), (double)(number_of_gridpoints+1));
            }
            else {
                w_upper = std::max(round2ffreq(wmax_in), (double)(number_of_gridpoints+1));
            }

#endif
        }
        w_upper = wmax_in;
        if (purely_positive) {
            w_lower = 0;
        }
        else {
            w_lower = -wmax_in;
        }

        if constexpr(!KELDYSH && !ZERO_T)
        {
#ifndef DENSEGRID
            // for Matsubara T>0: pick grid such that no frequencies occur twice
            if (type == 'b') {
                W_scale = wscale_from_wmax(W_scale, 2, w_upper, (number_of_gridpoints - 1) / 2);
            } else {
                W_scale = wscale_from_wmax(W_scale, 1, w_upper, number_of_gridpoints - 1);
            }
#else
            // for Matsubara T>0: pick grid such that no frequencies occur twice

        W_scale = Wscale_in;

#endif
        }
        else {
            W_scale = Wscale_in;
        }
    }
    void set_w_upper(freqType wmax) {double scale = W_scale; set_essential_parameters(wmax, scale);};
    void set_w_center(freqType w_center_in) {
        w_center = w_center_in;
        const double w_upper_centered = (w_upper - w_lower) * 0.5;
        w_upper = w_upper_centered + w_center_in;
        w_lower =-w_upper_centered + w_center_in;
        initialize_grid();
    }
    auto wscale_from_wmax(freqType & Wscale, freqType w1, freqType wmax, int N) -> double; // Not necessary if dense grid is chosen for Matsubara T>0

    void initialize_grid();
    void update_Wscale(freqType Wscale);

    /// core grid functionality (has to be super efficient)
    auto get_grid_index(freqType w_in) const -> int;
    int get_grid_index(freqType &t, freqType w_in) const;

    /// grid functions:
    auto t_from_frequency(freqType w) const -> freqType ;    // t(w)
    auto frequency_from_t(freqType t) const -> freqType;    // w(t)
};

template<>
class FrequencyGrid<hybridGrid> {
    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    friend class DataContainer;

public:
    /// essential grid parameters:
    int number_of_gridpoints;                       // total number of gridpoints
    freqType w_upper;                                 // largest positive frequency
    freqType w_lower;                                 // largest positive frequency
    freqType w_center = 0;        // symmetry point of frequency grid
    std::array<freqType,2> pos_section_boundaries;    // defines the sections [0, pos_section_boundaries[0]], [pos_section_boundaries[0], pos_section_boundaries[1]]... , [pos_section_boundaries[-1], w_upper]

    /// guess essential parameters from value of Lambda
    void guess_essential_parameters(double Lambda, const fRG_config& config);

//private:

    /// grid identifier:
    char type;
    int diag_class;
    bool purely_positive;

    /// auxiliary grid parameters:
    freqType t_upper;                    // largest point on auxiliary grid
    freqType t_lower;                    // smallest point on auxiliary grid
    std::array<double,2> aux_pos_section_boundaries;
    freqType spacing_auxiliary_gridpoint=1.;  // linear spacing on auxiliary grid for t
    freqType recip_curvature_quad;            // defines quadratic function f(t) = t^2 / recip_curvature_quad
    freqType recip_slope_lin;                 // defines linear function f(t) = pos_section_boundaries[0] + (t - aux_pos_section_boundaries[0]) / recip_slope_lin
    freqType factor_rat;                      // defines rational function f(t) = factor_rat / (1 - rescale_rat * t)
    freqType rescale_rat;                     // defines rational function f(t) = factor_rat / (1 - t / rescale_rat)
    void derive_auxiliary_parameters();     // derive auxiliary parameters from

    /// list of all frequencies:
    vec<freqType> all_frequencies;           // contains all frequencies w
    vec<freqType> auxiliary_grid;            // contains all t such that w(t) is the frequency

    /// grid functions:
    freqType frequency_from_t(freqType t) const;
    freqType t_from_frequency(freqType w) const;


public:
    /// constructor:
    FrequencyGrid(const char type_in, const unsigned int diag_class_in, const double Lambda, const fRG_config& config, const bool purely_positive=false) : type(type_in), diag_class(diag_class_in), purely_positive(purely_positive) {
#ifndef DENSEGRID
        assert(KELDYSH || ZERO_T);  // for Matsubara T>0 only allow dense grid
#endif
        guess_essential_parameters(Lambda, config);
        initialize_grid();
    }

    ///getter functions:
    freqType get_frequency(const int index) const {assert(index>=0); return all_frequencies[index];};
    freqType get_auxiliary_gridpoint(const int index) const {assert(index>=0); assert(index<number_of_gridpoints); return auxiliary_grid[index];};
    const vec<freqType>& get_all_frequencies() const {return all_frequencies;};
    const vec<freqType>& get_all_auxiliary_gridpoints() const {return auxiliary_grid;};
    freqType get_spacing_auxiliary_gridpoints() const {return spacing_auxiliary_gridpoint;}
    char get_type() const {return type;};
    char get_diag_class() const {return diag_class;};

    /// setter functions:
    void set_essential_parameters(freqType wmax_in, std::array<freqType,2> new_pos_section_boundaries) {
        w_upper = wmax_in;
        w_lower =-w_upper;
        pos_section_boundaries = std::move(new_pos_section_boundaries);
    }
    void set_w_upper(freqType wmax) {
        set_essential_parameters(wmax, pos_section_boundaries);
    }
    void initialize_grid();
    void update_pos_section_boundaries(std::array<freqType,2> new_pos_section_boundaries) {
        //pos_section_boundaries = std::move(new_pos_section_boundaries);
        pos_section_boundaries[0] = std::min(w_upper, new_pos_section_boundaries[0]);
        pos_section_boundaries[1] = std::max(pos_section_boundaries[0], new_pos_section_boundaries[1]);
        initialize_grid();
    }

    /// core grid functionality (has to be super efficient)
    int get_grid_index(freqType frequency)const;
    int get_grid_index(freqType& t, freqType frequency) const;


};


template<>
class FrequencyGrid<angularGrid> {
    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    friend class DataContainer;

public:
    /// essential grid parameters:
    int number_of_gridpoints;                       // total number of gridpoints
    double w_upper;                                 // largest  angle
    double w_lower;                                 // smallest angle
    freqType w_center = 0;        // symmetry point of frequency grid
    int number_of_intervals;    // defines the number of intervals between w_lower and w_upper on which we have quadratic functions
    double lin_fac = 1.e-4;      // "ratio" between linear and quadratic contribution in the grid function
    double power;
    double recip_power;

    /// guess essential parameters from value of Lambda
    void guess_essential_parameters(double Lambda, const fRG_config& config);

//private:

    /// grid identifier:
    char type;
    int diag_class;
    bool purely_positive;

    /// auxiliary grid parameters:
    double t_upper;                    // largest point on auxiliary grid
    double t_lower;                    // smallest point on auxiliary grid

    double spacing_auxiliary_gridpoint;  // linear spacing on auxiliary grid for t
    double half_of_interval_length_for_t;
    double half_of_interval_length_for_w;
    double half_of_interval_length_for_w_recip;
    double interval_length_for_w_recip;
    double quad_fac_recip;
    double lin_fac_to_power;
    void derive_auxiliary_parameters();     // derive auxiliary parameters from

    /// list of all frequencies:
    rvec all_frequencies;           // contains all frequencies w
    rvec auxiliary_grid;            // contains all t such that w(t) is the frequency

    /// grid functions:
    double frequency_from_t(double t) const;
    double t_from_frequency(double w) const;


public:
    /// constructor:
    FrequencyGrid(const char type_in, const unsigned int diag_class_in, const double Lambda, const fRG_config& config, const bool purely_positive) : type(type_in), diag_class(diag_class_in), purely_positive(purely_positive) {
#ifndef DENSEGRID
        assert(KELDYSH || ZERO_T);  // for Matsubara T>0 only allow dense grid
#else
        assert(false);
#endif
        guess_essential_parameters(Lambda, config);
        initialize_grid();
    }

    ///getter functions:
    double get_frequency(const int index) const {assert(index>=0); return all_frequencies[index];};
    double get_auxiliary_gridpoint(const int index) const {assert(index>=0); assert(index<number_of_gridpoints); return auxiliary_grid[index];};
    const rvec& get_all_frequencies() const {return all_frequencies;};
    const rvec& get_all_auxiliary_gridpoints() const {return auxiliary_grid;};
    double get_spacing_auxiliary_gridpoints() const {return spacing_auxiliary_gridpoint;}
    char get_type() const {return type;};
    char get_diag_class() const {return diag_class;};

    /// setter functions:
    void set_essential_parameters(double wmax_in, double number_of_intervals_in) {
        w_upper = wmax_in;
        w_lower = 0.;
        number_of_intervals = number_of_intervals_in;
    }
    void set_w_upper(double wmax) {
        set_essential_parameters(wmax, number_of_intervals);
    }
    void initialize_grid();
    void update_number_of_intervals(double number_of_intervals_in) {
        number_of_intervals = number_of_intervals_in;
        initialize_grid();
    }

    /// core grid functionality (has to be super efficient)
    int get_grid_index(double frequency)const;
    int get_grid_index(double& t, double frequency) const;

    void update_power(double power_in) {
        power = std::max(power_in, 1.); // don't go lower than 1
        initialize_grid();
    }

};

template<K_class k>
class bufferFrequencyGrid {
public:
#if GRID == 0
    using grid_type1 = FrequencyGrid<eliasGrid>;
    using grid_type2 = FrequencyGrid<eliasGrid>;
    using grid_type3 = FrequencyGrid<eliasGrid>;
#elif GRID == 1
     using grid_type1 = FrequencyGrid<hybridGrid>;
     using grid_type2 = FrequencyGrid<hybridGrid>;
     using grid_type3 = FrequencyGrid<hybridGrid>;
#else // GRID == 2
    using grid_type1 = FrequencyGrid<hybridGrid>;
    using grid_type2 = FrequencyGrid<k==k2||k==k2b||k==k3 ? angularGrid : hybridGrid>;
    using grid_type3 = FrequencyGrid<k==k2||k==k2b||k==k3 ? angularGrid : hybridGrid>;
#endif
    grid_type1   primary_grid;
    grid_type2 secondary_grid;
    grid_type3  tertiary_grid;
    double T;

    int get_diagclass() {
        if constexpr(k == selfenergy or k == k1) return 1;
        else if constexpr(k == k2 or k == k2b) return 2;
        else return 3;
    }

    bufferFrequencyGrid() :    primary_grid(k == selfenergy ? 'f' : 'b', get_diagclass(), 0, fRG_config(), GRID==2), secondary_grid('f', get_diagclass(), 0, fRG_config(), false), tertiary_grid('f', get_diagclass(), 0, fRG_config(), GRID==2) {};
    bufferFrequencyGrid(double Lambda, const fRG_config& config) :  T(config.T), primary_grid(k == selfenergy ? 'f' : 'b', get_diagclass(), Lambda, config, GRID==2 and (k==k2 or k==k2b or k==k3))
                                         , secondary_grid('f', get_diagclass(), Lambda, config, false)
                                         , tertiary_grid('f', get_diagclass(), Lambda, config, GRID==2)
    {};

    void guess_essential_parameters(double Lambda, const fRG_config& config) {
          primary_grid.guess_essential_parameters(Lambda, config);
        if constexpr(k != k1 and k != selfenergy) {
            secondary_grid.guess_essential_parameters(Lambda, config);
             tertiary_grid.guess_essential_parameters(Lambda, config);
        }
    }

    /// getter functions:
    auto get_freqGrid_b() const -> const grid_type1& {return   primary_grid;};
    auto get_freqGrid_f() const -> const grid_type2& {if constexpr(k != k1 and k != selfenergy)return secondary_grid; else assert(false);};//, "Exists no fermionic grid");};
    auto get_freqGrid_3() const -> const grid_type2& {if constexpr(k == k3 )return tertiary_grid; else assert(false);};//, "Exists no fermionic grid");};
    const freqType& get_wlower_b() const {return   primary_grid.w_lower;};
    const freqType& get_wupper_b() const {return   primary_grid.w_upper;};
    const freqType& get_wlower_f() const {if constexpr(k != k1 and k != selfenergy) return secondary_grid.w_lower; else assert(false);};//, "Exists no second grid");};
    const freqType& get_wupper_f() const {if constexpr(k != k1 and k != selfenergy) return secondary_grid.w_upper; else assert(false);};//, "Exists no second grid");};
    const freqType& get_tlower_b_aux() const {return   primary_grid.t_lower;};
    const freqType& get_tupper_b_aux() const {return   primary_grid.t_upper;};
    const freqType& get_tlower_f_aux() const {if constexpr(k != k1 and k != selfenergy) return secondary_grid.t_lower; else assert(false);};//, "Exists no second grid");};
    const freqType& get_tupper_f_aux() const {if constexpr(k != k1 and k != selfenergy) return secondary_grid.t_upper; else assert(false);};//, "Exists no second grid");};

    /// Check whether frequencies are in the frequency box:
    bool is_in_box(std::array<freqType,1> freqs) const {
        if constexpr(k == k1 or k == selfenergy) return std::abs(freqs[0] - primary_grid.w_center) <   primary_grid.w_upper + inter_tol;
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    bool is_in_box(std::array<freqType,2> freqs) const {
        K2_convert2internalFreqs(freqs[0], freqs[1]);
        if constexpr(k == k2) {
            if constexpr (!KELDYSH_FORMALISM and !ZERO_T) {
                return std::abs(freqs[0]) <   primary_grid.w_upper + inter_tol and (freqs[1] > -secondary_grid.w_upper - inter_tol and freqs[1] < secondary_grid.w_upper + signFlipCorrection_MF(freqs[0]) + inter_tol );

            }
            else {
                return std::abs(freqs[0] - primary_grid.w_center) <   primary_grid.w_upper + inter_tol and std::abs(freqs[1] - secondary_grid.w_center) < secondary_grid.w_upper + inter_tol;
            }

        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    bool is_in_box(std::array<freqType,3> freqs) const {
        K3_convert2internalFreqs(freqs[0], freqs[1], freqs[2]);
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION)) return std::abs(freqs[0] - primary_grid.w_center) <   primary_grid.w_upper + inter_tol and std::abs(freqs[1] - secondary_grid.w_center) < secondary_grid.w_upper + inter_tol and std::abs(freqs[2] - tertiary_grid.w_center) < tertiary_grid.w_upper + inter_tol;
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    /// determine the grid indices for the frequencies + determine normalized distance to next smaller grid point on all_frequencies
    void get_grid_index(std::array<my_index_t,1>& idx, std::array<double,1>& dw_normalized, const std::array<freqType,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            freqType w = freqs[0];
            int iw =   primary_grid.get_grid_index(w);
            idx[0] = iw;
            double w_low =   primary_grid.get_frequency(iw);
            double w_high=   primary_grid.get_frequency(iw+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_grid_index(std::array<my_index_t,2>& idx, std::array<double,2>& dw_normalized, const std::array<freqType,2>& freqs) const {
        if constexpr(k == k2)  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            K2_convert2internalFreqs(w, v);

            int iw =   primary_grid.get_grid_index(w);
            int iv = secondary_grid.get_grid_index(v);
            idx[0] = iw;
            idx[1] = iv;
            double w_low =   primary_grid.get_frequency(iw);
            double w_high=   primary_grid.get_frequency(iw+1);
            double v_low = secondary_grid.get_frequency(iv);
            double v_high= secondary_grid.get_frequency(iv+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
            dw_normalized[1] = (v - v_low) / (v_high - v_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_grid_index(std::array<my_index_t,3>& idx, std::array<double,3>& dw_normalized, const std::array<freqType,3>& freqs) const {
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION))  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            freqType vp = freqs[2];
            K3_convert2internalFreqs(w, v, vp);

            const int iw =   primary_grid.get_grid_index(freqs[0]);
            const int iv = secondary_grid.get_grid_index(freqs[1]);
            const int ivp=  tertiary_grid.get_grid_index(freqs[2]);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            const double w_low =   primary_grid.get_frequency(iw);
            const double w_high=   primary_grid.get_frequency(iw+1);
            const double v_low = secondary_grid.get_frequency(iv);
            const double v_high= secondary_grid.get_frequency(iv+1);
            const double vp_low =tertiary_grid.get_frequency(ivp);
            const double vp_high=tertiary_grid.get_frequency(ivp+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
            dw_normalized[1] = (v - v_low) / (v_high - v_low);
            dw_normalized[2] = (vp-vp_low) / (vp_high-vp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    /// determine the grid indices for the frequencies + determine normalized distance to next smaller grid point on auxiliary grid
    void get_auxgrid_index(std::array<my_index_t,1>& idx, std::array<double,1>& dt_normalized, const std::array<freqType,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            //double w = freqs[0];
            freqType tw;
            int iw =   primary_grid.get_grid_index(tw, freqs[0]);
            idx[0] = iw;
            double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_auxgrid_index(std::array<my_index_t,2>& idx, std::array<freqType,2>& dt_normalized, const std::array<freqType,2>& freqs) const {
        if constexpr(k == k2)  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            K2_convert2internalFreqs(w, v);

            freqType tw, tv;
            int iw =   primary_grid.get_grid_index(tw, w);
            int iv = secondary_grid.get_grid_index(tv, v);
            idx[0] = iw;
            idx[1] = iv;
            double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            double tv_low = secondary_grid.get_auxiliary_gridpoint(iv);
            double tv_high= secondary_grid.get_auxiliary_gridpoint(iv+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
            dt_normalized[1] = (tv - tv_low) / (tv_high - tv_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_auxgrid_index(std::array<my_index_t,3>& idx, std::array<double,3>& dt_normalized, const std::array<freqType,3>& freqs) const {
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION))  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            freqType vp = freqs[2];
            K3_convert2internalFreqs(w, v, vp);

            freqType tw, tv, tvp;
            int iw =   primary_grid.get_grid_index(tw, w);
            int iv = secondary_grid.get_grid_index(tv, v);
            int ivp=  tertiary_grid.get_grid_index(tvp,vp);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            double tv_low = secondary_grid.get_auxiliary_gridpoint(iv);
            double tv_high= secondary_grid.get_auxiliary_gridpoint(iv+1);
            double tvp_low =tertiary_grid.get_auxiliary_gridpoint(ivp);
            double tvp_high=tertiary_grid.get_auxiliary_gridpoint(ivp+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
            dt_normalized[1] = (tv - tv_low) / (tv_high - tv_low);
            dt_normalized[2] = (tvp-tvp_low) / (tvp_high-tvp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    /// determine the grid indices for the frequencies + determine UN-normalized distance to next smaller grid point on auxiliary grid
    void get_auxgrid_index_unnormalized(std::array<my_index_t,1>& idx, std::array<double,1>& dt_unnormalized, const std::array<freqType,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            //double w = freqs[0];
            freqType tw;
            int iw =   primary_grid.get_grid_index(tw, freqs[0]);
            idx[0] = iw;
            double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_auxgrid_index_unnormalized(std::array<my_index_t,2>& idx, std::array<double,2>& dt_unnormalized, const std::array<freqType,2>& freqs) const {
        if constexpr(k == k2)  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            K2_convert2internalFreqs(w, v);

            freqType tw, tv;
            int iw =   primary_grid.get_grid_index(tw, w);
            int iv = secondary_grid.get_grid_index(tv, v);
            idx[0] = iw;
            idx[1] = iv;
            double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            //double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            double tv_low = secondary_grid.get_auxiliary_gridpoint(iv);
            //double tv_high= secondary_grid.get_auxiliary_gridpoint(iv+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
            dt_unnormalized[1] = (tv - tv_low); // / (tv_high - tv_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_auxgrid_index_unnormalized(std::array<my_index_t,3>& idx, std::array<double,3>& dt_unnormalized, const std::array<freqType,3>& freqs) const {
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION))  {
            freqType w  = freqs[0];
            freqType v  = freqs[1];
            freqType vp = freqs[2];
            K3_convert2internalFreqs(w, v, vp);

            freqType tw, tv, tvp;
            const int iw =   primary_grid.get_grid_index(tw, w);
            const int iv = secondary_grid.get_grid_index(tv, v);
            const int ivp=  tertiary_grid.get_grid_index(tvp,vp);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            const double tw_low =   primary_grid.get_auxiliary_gridpoint(iw);
            //const double tw_high=   primary_grid.get_auxiliary_gridpoint(iw+1);
            const double tv_low = secondary_grid.get_auxiliary_gridpoint(iv);
            //const double tv_high= secondary_grid.get_auxiliary_gridpoint(iv+1);
            const double tvp_low = tertiary_grid.get_auxiliary_gridpoint(ivp);
            //const double tvp_high= tertiary_grid.get_auxiliary_gridpoint(ivp+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
            dt_unnormalized[1] = (tv - tv_low); // / (tv_high - tv_low);
            dt_unnormalized[2] = (tvp-tvp_low); // / (tvp_high-tvp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    /// get frequencies corresponding to a set of grid indices
    void get_freqs_w(freqType &w, const int iw) const {
        if constexpr(k == k1 or k == selfenergy) w =   primary_grid.get_frequency(iw);
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_w(freqType &w, freqType &v, const int iw, const int iv) const {
        if constexpr(k == k2)
        {
            w =   primary_grid.get_frequency(iw);
            v = secondary_grid.get_frequency(iv);
            K2_convert2naturalFreqs(w, v);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_w(freqType &w, freqType &v,  freqType &vp, const int iw, const int iv, const int ivp) const {
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION)) {
            w =   primary_grid.get_frequency(iw);
            v = secondary_grid.get_frequency(iv);
            vp = tertiary_grid.get_frequency(ivp);
            K3_convert2naturalFreqs(w, v, vp);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    // wraps above functions for arrays as input
    template<size_t num>
    void get_freqs_w(std::array<freqType, num>& freqs, std::array<my_index_t , num>& i_freqs) const{
        if constexpr(num == 1) {
            get_freqs_w(freqs[0], i_freqs[0]);
        }
        else if constexpr(num == 2) {
            get_freqs_w(freqs[0], freqs[1], i_freqs[0], i_freqs[1]);
        }
        else if constexpr(num == 3) {
            get_freqs_w(freqs[0], freqs[1], freqs[2], i_freqs[0], i_freqs[1], i_freqs[2]);
        }
        else assert(false);
    }


    /// only used in test_functions:
    void get_freqs_aux(freqType &w, const int iw) const {
        if constexpr(k == k1 or k == selfenergy) w =   primary_grid.get_auxiliary_gridpoint(iw);
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_aux(freqType &w, freqType &v, const int iw, const int iv) const {
        if constexpr(k == k2) {
            w =   primary_grid.get_auxiliary_gridpoint(iw);
            v = secondary_grid.get_auxiliary_gridpoint(iv);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_aux(freqType &w, freqType &v, double &vp, const int iw, const int iv, const int ivp) const {
        if constexpr(k == k3 or (k == k2 and SBE_DECOMPOSITION)) {
            w =   primary_grid.get_auxiliary_gridpoint(iw);
            v = secondary_grid.get_auxiliary_gridpoint(iv);
            vp=  tertiary_grid.get_auxiliary_gridpoint(ivp);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    auto gridtransf_f(freqType w) const -> freqType {if constexpr(k != k1 and k != selfenergy) return secondary_grid.t_from_frequency(w); else assert(false);};//, "Exists no second grid");};

};

/*******************************************    FREQUENCY GRID    *****************************************************/

double grid_transf_lin(double w, double W_scale);
//int    grid_transf_lin(int w, int W_scale);
double grid_transf_v1(double w, double W_scale, double w_center);
double grid_transf_v2(double w, double W_scale, double w_center);
double grid_transf_v3(double w, double W_scale, double w_center);
double grid_transf_v4(double w, double W_scale, double w_center);
double grid_transf_log(double w, double W_scale);
freqType grid_transf_inv_lin(freqType W, freqType W_scale);
double grid_transf_inv_v1(double t, double W_scale, double w_center);
double grid_transf_inv_v2(double t, double W_scale, double w_center);
double grid_transf_inv_v3(double t, double W_scale, double w_center);
double grid_transf_inv_v4(double t, double W_scale, double w_center);
double grid_transf_inv_log(double t, double W_scale);
freqType wscale_from_wmax_v1(freqType & Wscale, freqType w1, freqType wmax, int N);
freqType wscale_from_wmax_v2(freqType & Wscale, freqType w1, freqType wmax, int N);
freqType wscale_from_wmax_v3(freqType & Wscale, freqType w1, freqType wmax, int N);
freqType wscale_from_wmax_lin(freqType & Wscale, freqType w1, freqType wmax, int N);



#endif //KELDYSH_MFRG_FREQUENCY_GRID_HPP
