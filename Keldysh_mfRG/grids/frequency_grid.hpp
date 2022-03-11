/**
 * Set up the frequency grid //
 * Functions that initialize the grid and provide conversion between doubles and grid indices.
 * Three different grid types:
 * GRID=1: log grid -- to be implemented
 * GRID=2: linear grid
 * GRID=3: non-linear grid in different versions:
 *     version 1: w(t) = W_scale t  /sqrt(1 - t^2)
 *     version 2: w(t) = W_scale t^2/sqrt(1 - t^2)
 *     version 3: w(t) = W_scale t  /    (1 - t^2)
 *     version 4: w(t) = W_scale t^2/    (1 - t^2)
 * The grid points are obtained as follows:
 *  The values of t are evenly distributed between t_lower and t_upper on a sub-interval of [-1,1].
 *  The frequencies w(t) are obtained according to the function.
 *  W_scale determines 'how non-linear' the grid behaves.
 * GRID=4: tangent grid w = c1 * tan( c2 * i + c3 )
 */

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_HPP
#define KELDYSH_MFRG_FREQUENCY_GRID_HPP

#include <cmath>        // for sqrt, log, exp
#include "../data_structures.hpp"
#include "../utilities/util.hpp"
#include "../parameters/master_parameters.hpp" // for frequency/Lambda limits and number of frequency/Lambda points
#include "../utilities/math_utils.hpp"
#include <cassert>
#include "H5Cpp.h"

// TODO(low): implement functions used for GRID=3 also for GRID=1,2,4
/// TODO: implement hybrid frequency grid, improve/unify treatment of different frequency meshes and their frequency parameters

template<K_class k, typename Q> class vertexDataContainer; // forward declaration
template<typename Q> class State; // forward declaration

#define PARAMETRIZED_GRID
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
#define DENSEGRID
#endif
#ifdef DENSEGRID
const bool dense = true;
#else
const bool dense = false;
#endif

class FrequencyGrid {
    template<K_class k, typename Q> friend class vertexDataContainer;
    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type> friend class DataContainer;
    //friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);
    friend void init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, FrequencyGrid& freqgrid, int Lambda_it);

    char type;
    unsigned int diag_class;
    rvec ws;                     // frequency grid
    rvec ts;                    // linear auxiliary grid (related to ws by ws(ts)=grid_transf_inv(ts))
public:
    int N_w;
    double w_upper;             // lower bound of frequency grid
    double w_lower;             // lower bound of frequency grid
    double t_upper;             // upper bound of auxiliary grid
    double t_lower;             // lower bound of auxiliary grid
    double W_scale;             // non-linearity of grid_transf()
    double dt;                  // spacing on linear auxiliary grid
    double U_factor = 0./3.;   // determines scale_factor()
    double Delta_factor = 10.;  // determines scale_factor()

    /**
     * This constructor initializes a frequency grid with the global values. This is not needed anymore!
     * @param type_in
     * @param diag_class_in
     * @param Lambda
     */
    FrequencyGrid(char type_in, unsigned int diag_class_in, double Lambda) : type(type_in), diag_class(diag_class_in) {
        switch (type) {
            case 'b':
                switch (diag_class) {
                    case 1:
                        N_w = nBOS;
                        if (KELDYSH) {
                            U_factor = 0. / 3.;
                            Delta_factor = 5.;
                        }
                        else {
                            U_factor = 40./3.;
                            Delta_factor = 40.;
                        }
                        break;
                    case 2:
                        N_w = nBOS2;
#ifdef ROTATEK2
                        if (KELDYSH){
                            U_factor = 10./3.;
                            Delta_factor = 10.;
                        }
                        else{
                            U_factor = 10./3.;
                            Delta_factor = 10.;
                        }
#else
                        if (KELDYSH){
                            U_factor = 0./3.;
                            Delta_factor = 15.;
                        }
                        else{
                            U_factor = 10./3.;
                            Delta_factor = 10.;
                        }
#endif
                        break;
                    case 3:
                        N_w = nBOS3;
                        break;
                    default:;
                }
                break;
            case 'f':
                switch (diag_class) {
                    case 1:
                        N_w = nFER;
                        if (KELDYSH) {
                            U_factor = 0. / 3.;
                            Delta_factor = 10.;
                        }
                        else {
                            U_factor = 2./3.;
                            Delta_factor = 2.;
                        }
                        if (HUBBARD_MODEL){ //TODO(medium): Just a hotfix for the Hubbard model. Avoids that one runs out of the frequency box when integrating for the bubble.
                            U_factor *= 1.5;
                            Delta_factor *= 1.5;
                        }
                        break;
                    case 2:
                        N_w = nFER2;
#ifdef ROTATEK2
                        /// Needs to be the same as for 'b'!!!
                        if (KELDYSH) {
                            U_factor = 10. / 3.;
                            Delta_factor = 10.;
                        }
                        else {
                            U_factor = 10./3.;
                            Delta_factor = 10.;
                        }
#else
                        if (KELDYSH) {
                            U_factor = 0. / 3.;
                            Delta_factor = 20.;
                        }
                        else {
                            U_factor = 4./3.;
                            Delta_factor = 4.;
                        }
#endif
                        break;
                    case 3:
                        N_w = nFER3;
                        break;
                    default:;
                }
                break;
            default:;
        }
        ws = rvec (N_w);
        ts = rvec (N_w);

        rescale_grid(Lambda);
    };



    auto operator= (const FrequencyGrid& freqGrid) -> FrequencyGrid& {
        assert(this->type == freqGrid.type);
        assert(this->diag_class == freqGrid.diag_class);
        this->N_w = freqGrid.N_w;
        this->w_upper = freqGrid.w_upper;
        this->w_lower = freqGrid.w_lower;
        this->t_upper = freqGrid.t_upper;
        this->t_lower = freqGrid.t_lower;
        this->W_scale = freqGrid.W_scale;
        this->dt = freqGrid.dt;
        this->U_factor = freqGrid.U_factor;
        this->Delta_factor = freqGrid.Delta_factor;
        this->ws = freqGrid.ws;
        this->ts = freqGrid.ts;
        return *this;
    }

    int get_diag_class() const {return diag_class;}
    char get_type() const {return type;}
    auto get_ws(int index) const -> double {assert(index>=0); assert(index<N_w); assert(isfinite(ws[index])); return ws[index];};
    auto get_ts(int index) const -> double {assert(index>=0); assert(index<N_w); assert(isfinite(ts[index])); return ts[index];};
    auto get_ws_vec() const -> vec<double> {return ws;}
    auto get_ts_vec() const -> vec<double> {return ts;}
    auto scale_factor(double Lambda) -> double;
    void initialize_grid();
    void set_W_scale(double scale);
    void set_w_upper(double wmax);
    void rescale_grid(double Lambda);
    void update_Wscale(double Wscale);
    auto fconv(double w_in, bool safety=dense) const -> int;
    auto grid_transf(double w) const -> double;
    auto grid_transf_inv(double t) const -> double;
    auto wscale_from_wmax(double & Wscale, double w1, double wmax, int N) -> double;

    int fconv(double &t, double w_in) const;
};


enum freqGrid_identifier {grid4selfenergy, grid4K1, grid4K2, grid4K3};

/**
 * Initializes frequency grids for a vertex
 */
template<K_class>
class VertexFrequencyGrid {};
template<>
class VertexFrequencyGrid<k1> {
public:
    FrequencyGrid b;

    VertexFrequencyGrid<k1>() :  b('b', 1, 0) {};
    VertexFrequencyGrid<k1>(double Lambda) : b('b', 1, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
    }

    void get_freqs_w(double &w, const int iw) const {
        w = b.get_ws(iw);
    }

    void get_freqs_aux(double &w, const int iw) const {
        w = b.get_ts(iw);
    }

    auto get_freqGrid_b() const -> FrequencyGrid {return b;};
//
    const double& get_wlower_b() const {return b.w_lower;};
    const double& get_wupper_b() const {return b.w_upper;};
    const double& get_tlower_aux() const {return b.t_lower;};
    const double& get_tupper_aux() const {return b.t_upper;};

    // currently only used in test_interpolations:
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
//
    //void get_freqs_w_b(double& w, int i) const {w = b.ws[i];};
    //void get_freqs_aux_b(double& w, int iw) const {w = b.ts[iw];};
};
template<>
class VertexFrequencyGrid<k2> {
public:
    FrequencyGrid b;
    FrequencyGrid f;


    VertexFrequencyGrid<k2>() :  b('b', 2, 0), f('f', 2, 0) {};
    VertexFrequencyGrid<k2>(double Lambda) : b('b', 2, Lambda),
                                             f('f', 2, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
        f.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
        f.set_W_scale(scale);
        f.set_w_upper(scale*15.);
        f.initialize_grid();
    }

    auto get_freqGrid_b() const -> const FrequencyGrid& {return b;};
    auto get_freqGrid_f() const -> const FrequencyGrid& {return f;};
//
    const double& get_wlower_b() const {return b.w_lower;};
    const double& get_wupper_b() const {return b.w_upper;};
    const double& get_wlower_f() const {return f.w_lower;};
    const double& get_wupper_f() const {return f.w_upper;};
    const double& get_tlower_b_aux() const {return b.t_lower;};
    const double& get_tupper_b_aux() const {return b.t_upper;};
    const double& get_tlower_f_aux() const {return f.t_lower;};
    const double& get_tupper_f_aux() const {return f.t_upper;};
    auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    auto gridtransf_f(double w) const -> double {return f.grid_transf(w);};
    auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    auto gridtransf_inv_f(double t) const -> double {return f.grid_transf_inv(t);};
//
    void get_freqs_w(double &w, double &v, const int iw, const int iv) const {
        w = b.get_ws(iw);
        v = f.get_ws(iv);
        K2_convert2naturalFreqs(w, v);

    }

    void get_freqs_aux(double &w, double &v, const int iw, const int iv) const {
        w = b.get_ts(iw);
        v = f.get_ts(iv);
    }
};
template<>
class VertexFrequencyGrid<k3> {
public:
    FrequencyGrid b;
    FrequencyGrid f;


    VertexFrequencyGrid<k3>() : b('b', 3, 0), f('f', 3, 0) {};
    VertexFrequencyGrid<k3>(double Lambda) : b('b', 3, Lambda),
                                             f('f', 3, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
        f.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
        f.set_W_scale(scale);
        f.set_w_upper(scale*15.);
        f.initialize_grid();
    }

    auto get_freqGrid_b() const -> FrequencyGrid {return b;};
    auto get_freqGrid_f() const -> FrequencyGrid {return f;};

    const double& get_wlower_b() const {return b.w_lower;};
    const double& get_wupper_b() const {return b.w_upper;};
    const double& get_wlower_f() const {return f.w_lower;};
    const double& get_wupper_f() const {return f.w_upper;};
    const double& get_tlower_b_aux() const {return b.t_lower;};
    const double& get_tupper_b_aux() const {return b.t_upper;};
    const double& get_tlower_f_aux() const {return f.t_lower;};
    const double& get_tupper_f_aux() const {return f.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_f(double w) const -> double {return f.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    //auto gridtransf_inv_f(double t) const -> double {return f.grid_transf_inv(t);};
//
    void get_freqs_w(double &w, double &v, double& vp, const int iw, const int iv, const int ivp, const char channel) const {
        w = b.get_ws(iw);
        v = f.get_ws(iv);
        vp= f.get_ws(ivp);

        if (BOSONIC_PARAM_FOR_K3) {
            if (channel == 'a') { switch2naturalFreqs<'a'>(w, v, vp); }
            else if (channel == 'p') { switch2naturalFreqs<'p'>(w, v, vp); }
            else if (channel == 't') { switch2naturalFreqs<'t'>(w, v, vp); }
        }
    }

    void get_freqs_aux(double &w, double &v, double& vp, const int iw, const int iv, const int ivp) const {
        w = b.get_ts(iw);
        v = f.get_ts(iv);
        vp= f.get_ts(ivp);
    }
};



template<K_class k>
class bufferFrequencyGrid {
public:
    FrequencyGrid b;
    FrequencyGrid f;

    int get_diagclass() {
        if constexpr(k == selfenergy or k == k1) return 1;
        else if constexpr(k == k2 or k == k2b) return 2;
        else return 3;
    }

    bufferFrequencyGrid() :  b(k == selfenergy ? 'f' : 'b', get_diagclass(), 0), f('f', get_diagclass(), 0) {};
    bufferFrequencyGrid(double Lambda) : b(k == selfenergy ? 'f' : 'b', get_diagclass(), Lambda), f('f', get_diagclass(), Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
        if constexpr(k != k1 and k != selfenergy)f.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {
            b.set_W_scale(scale);
            b.set_w_upper(scale*15.);
            b.initialize_grid();
        if constexpr(k != k1 and k != selfenergy) {
            f.set_W_scale(scale);
            f.set_w_upper(scale*15.);
            f.initialize_grid();
        }
    }

    auto get_freqGrid_b() const -> const FrequencyGrid& {return b;};
    auto get_freqGrid_f() const -> const FrequencyGrid& {if constexpr(k != k1 and k != selfenergy)return f; else assert(false);};//, "Exists no fermionic grid");};
//
    const double& get_wlower_b() const {return b.w_lower;};
    const double& get_wupper_b() const {return b.w_upper;};
    const double& get_wlower_f() const {if constexpr(k != k1 and k != selfenergy) return f.w_lower; else assert(false);};//, "Exists no second grid");};
    const double& get_wupper_f() const {if constexpr(k != k1 and k != selfenergy) return f.w_upper; else assert(false);};//, "Exists no second grid");};
    const double& get_tlower_b_aux() const {return b.t_lower;};
    const double& get_tupper_b_aux() const {return b.t_upper;};
    const double& get_tlower_f_aux() const {if constexpr(k != k1 and k != selfenergy) return f.t_lower; else assert(false);};//, "Exists no second grid");};
    const double& get_tupper_f_aux() const {if constexpr(k != k1 and k != selfenergy) return f.t_upper; else assert(false);};//, "Exists no second grid");};
    auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    auto gridtransf_f(double w) const -> double {if constexpr(k != k1 and k != selfenergy) return f.grid_transf(w); else assert(false);};//, "Exists no second grid");};
    auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    auto gridtransf_inv_f(double t) const -> double {if constexpr(k != k1 and k != selfenergy) return f.grid_transf_inv(t); else assert(false);};//, "Exists no second grid");};
//
    bool is_in_box(std::array<double,1> freqs) const {
        if constexpr(k == k1 or k == selfenergy) return std::abs(freqs[0]) < b.w_upper + inter_tol;
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    bool is_in_box(std::array<double,2> freqs) const {
        if constexpr(k == k2) return std::abs(freqs[0]) < b.w_upper + inter_tol and std::abs(freqs[1]) < f.w_upper + inter_tol;
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    bool is_in_box(std::array<double,3> freqs) const {
        if constexpr(k == k3) return std::abs(freqs[0]) < b.w_upper + inter_tol and std::abs(freqs[1]) < f.w_upper + inter_tol and std::abs(freqs[2]) < f.w_upper + inter_tol;
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    void fconv(std::array<my_index_t,1>& idx, std::array<double,1>& dw_normalized, const std::array<double,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            double w = freqs[0];
            int iw = b.fconv(w);
            idx[0] = iw;
            double w_low = b.get_ws(iw);
            double w_high= b.get_ws(iw+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv(std::array<my_index_t,2>& idx, std::array<double,2>& dw_normalized, const std::array<double,2>& freqs) const {
        if constexpr(k == k2)  {
            double w  = freqs[0];
            double v  = freqs[1];
            int iw = b.fconv(freqs[0]);
            int iv = f.fconv(freqs[1]);
            idx[0] = iw;
            idx[1] = iv;
            double w_low = b.get_ws(iw);
            double w_high= b.get_ws(iw+1);
            double v_low = f.get_ws(iv);
            double v_high= f.get_ws(iv+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
            dw_normalized[1] = (v - v_low) / (v_high - v_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv(std::array<my_index_t,3>& idx, std::array<double,3>& dw_normalized, const std::array<double,3>& freqs) const {
        if constexpr(k == k3)  {
            double w  = freqs[0];
            double v  = freqs[1];
            double vp = freqs[2];
            int iw = b.fconv(freqs[0]);
            int iv = f.fconv(freqs[1]);
            int ivp= f.fconv(freqs[2]);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            double w_low = b.get_ws(iw);
            double w_high= b.get_ws(iw+1);
            double v_low = f.get_ws(iv);
            double v_high= f.get_ws(iv+1);
            double vp_low =f.get_ws(ivp);
            double vp_high=f.get_ws(ivp+1);
            dw_normalized[0] = (w - w_low) / (w_high - w_low);
            dw_normalized[1] = (v - v_low) / (v_high - v_low);
            dw_normalized[2] = (vp-vp_low) / (vp_high-vp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    void fconv_on_aux(std::array<my_index_t,1>& idx, std::array<double,1>& dt_normalized, const std::array<double,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            //double w = freqs[0];
            double tw;
            int iw = b.fconv(tw, freqs[0]);
            idx[0] = iw;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv_on_aux(std::array<my_index_t,2>& idx, std::array<double,2>& dt_normalized, const std::array<double,2>& freqs) const {
        if constexpr(k == k2)  {
            //double w  = freqs[0];
            //double v  = freqs[1];
            double tw, tv;
            int iw = b.fconv(tw, freqs[0]);
            int iv = f.fconv(tv, freqs[1]);
            idx[0] = iw;
            idx[1] = iv;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            double tv_low = f.get_ts(iv);
            double tv_high= f.get_ts(iv+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
            dt_normalized[1] = (tv - tv_low) / (tv_high - tv_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv_on_aux(std::array<my_index_t,3>& idx, std::array<double,3>& dt_normalized, const std::array<double,3>& freqs) const {
        if constexpr(k == k3)  {
            //double w  = freqs[0];
            //double v  = freqs[1];
            //double vp = freqs[2];
            double tw, tv, tvp;
            int iw = b.fconv(tw, freqs[0]);
            int iv = f.fconv(tv, freqs[1]);
            int ivp= f.fconv(tvp,freqs[2]);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            double tv_low = f.get_ts(iv);
            double tv_high= f.get_ts(iv+1);
            double tvp_low =f.get_ts(ivp);
            double tvp_high=f.get_ts(ivp+1);
            dt_normalized[0] = (tw - tw_low) / (tw_high - tw_low);
            dt_normalized[1] = (tv - tv_low) / (tv_high - tv_low);
            dt_normalized[2] = (tvp-tvp_low) / (tvp_high-tvp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv_on_aux_unnormalized(std::array<my_index_t,1>& idx, std::array<double,1>& dt_unnormalized, const std::array<double,1>& freqs) const {
        if constexpr(k == k1 or k == selfenergy)  {
            //double w = freqs[0];
            double tw;
            int iw = b.fconv(tw, freqs[0]);
            idx[0] = iw;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv_on_aux_unnormalized(std::array<my_index_t,2>& idx, std::array<double,2>& dt_unnormalized, const std::array<double,2>& freqs) const {
        if constexpr(k == k2)  {
            //double w  = freqs[0];
            //double v  = freqs[1];
            double tw, tv;
            int iw = b.fconv(tw, freqs[0]);
            int iv = f.fconv(tv, freqs[1]);
            idx[0] = iw;
            idx[1] = iv;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            double tv_low = f.get_ts(iv);
            double tv_high= f.get_ts(iv+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
            dt_unnormalized[1] = (tv - tv_low); // / (tv_high - tv_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void fconv_on_aux_unnormalized(std::array<my_index_t,3>& idx, std::array<double,3>& dt_unnormalized, const std::array<double,3>& freqs) const {
        if constexpr(k == k3)  {
            //double w  = freqs[0];
            //double v  = freqs[1];
            //double vp = freqs[2];
            double tw, tv, tvp;
            int iw = b.fconv(tw, freqs[0]);
            int iv = f.fconv(tv, freqs[1]);
            int ivp= f.fconv(tvp,freqs[2]);
            idx[0] = iw;
            idx[1] = iv;
            idx[2] = ivp;
            double tw_low = b.get_ts(iw);
            double tw_high= b.get_ts(iw+1);
            double tv_low = f.get_ts(iv);
            double tv_high= f.get_ts(iv+1);
            double tvp_low =f.get_ts(ivp);
            double tvp_high=f.get_ts(ivp+1);
            dt_unnormalized[0] = (tw - tw_low); // / (tw_high - tw_low);
            dt_unnormalized[1] = (tv - tv_low); // / (tv_high - tv_low);
            dt_unnormalized[2] = (tvp-tvp_low); // / (tvp_high-tvp_low);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    void get_freqs_w(double &w, const int iw) const {
        if constexpr(k == k1 or k == selfenergy) w = b.get_ws(iw);
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_w(double &w, double &v, const int iw, const int iv) const {
        if constexpr(k == k2)
        {
            w = b.get_ws(iw);
            v = f.get_ws(iv);
            //K2_convert2naturalFreqs(w, v);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_w(double &w, double &v,  double &vp, const int iw, const int iv, const int ivp) const {
        if constexpr(k == k3) {
            w = b.get_ws(iw);
            v = f.get_ws(iv);
            vp = f.get_ws(ivp);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }

    template<size_t num>
    void get_freqs_w(std::array<double, num>& freqs, std::array<my_index_t , num>& i_freqs) const{
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

    void get_freqs_aux(double &w, const int iw) const {
        if constexpr(k == k1 or k == selfenergy) w = b.get_ts(iw);
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_aux(double &w, double &v, const int iw, const int iv) const {
        if constexpr(k == k2) {
            w = b.get_ts(iw);
            v = f.get_ts(iv);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
    void get_freqs_aux(double &w, double &v, double &vp, const int iw, const int iv, const int ivp) const {
        if constexpr(k == k3) {
            w = b.get_ts(iw);
            v = f.get_ts(iv);
            vp= f.get_ts(ivp);
        }
        else assert(false); // "Inconsistent number of frequency arguments.");
    }
};

/*******************************************    FREQUENCY GRID    *****************************************************/


#if GRID==1
/***********************************************    LOG GRID    *******************************************************/

double sgn(double x);
double grid_transf_b(double w);
double grid_transf_b_inv(double W);
double grid_transf_f(double w);
double grid_transf_f_inv(double W);
void setUpBosGrid(rvec& freqs, int nfreqs);
void setUpFerGrid(rvec& freqs, int nfreqs);
auto fconv_bos(double w, int nfreqs) -> int;
auto fconv_fer(double w, int nfreqs) -> int;


#elif GRID==2
/*********************************************    LINEAR GRID    ******************************************************/

void setUpBosGrid(rvec& freqs, int nfreqs);
void setUpFerGrid(rvec& freqs, int nfreqs);
auto fconv_bos(double w, int nfreqs) -> int;
auto fconv_fer(double v, int nfreqs) -> int;

#include <tuple>   // return several indices
auto fconv_K1_a(double w) -> int;
auto fconv_K2_a(double w, double v1) -> tuple<int, int>;
auto fconv_K3_a(double w, double v1, double v2) -> tuple<int, int, int>;
auto fconv_K1_p(double w) -> int;
auto fconv_K2_p(double w, double v1) -> tuple<int, int>;
auto fconv_K3_p(double w, double v1, double v2) -> tuple<int, int, int>;
auto fconv_K1_t(double w) -> int;
auto fconv_K2_t(double w, double v1) -> tuple<int, int>;
auto fconv_K3_t(double w, double v1, double v2) -> tuple<int, int, int>;


#elif GRID==3
/*******************************************    NON-LINEAR GRID    ****************************************************/

double sgn(double x);

double grid_transf_lin(double w, double W_scale);
double grid_transf_v1(double w, double W_scale);
double grid_transf_v2(double w, double W_scale);
double grid_transf_v3(double w, double W_scale);
double grid_transf_v4(double w, double W_scale);
double grid_transf_log(double w, double W_scale);
double grid_transf_inv_lin(double W, double W_scale);
double grid_transf_inv_v1(double t, double W_scale);
double grid_transf_inv_v2(double t, double W_scale);
double grid_transf_inv_v3(double t, double W_scale);
double grid_transf_inv_v4(double t, double W_scale);
double grid_transf_inv_log(double t, double W_scale);
double wscale_from_wmax_v1(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v2(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v3(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_lin(double & Wscale, double w1, double wmax, int N);
double integration_measure_v1(double t, double W_scale);
double integration_measure_v2(double t, double W_scale);
double integration_measure_v3(double t, double W_scale);
double integration_measure_v4(double t, double W_scale);

namespace freqGrid {
    auto shrink_freq_box(const FrequencyGrid& freqGrid, const double  rel_tail_threshold, const vec<double>& maxabs_along_x, bool verbose=true) -> FrequencyGrid;
}


#elif GRID==4
/*********************************************    TAN GRID    ******************************************************/

void setUpBosGrid(rvec& freqs, int nfreqs);
void setUpFerGrid(rvec& freqs, int nfreqs);
auto fconv_bos(double w, int nfreqs) -> int;
auto fconv_fer(const double v, int nfreqs) -> int;

#endif // GRID

#endif //KELDYSH_MFRG_FREQUENCY_GRID_HPP
