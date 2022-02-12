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
    double U_factor = 10./3.;   // determines scale_factor()
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
                            U_factor = 30. / 3.;
                            Delta_factor = 30.;
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
                            U_factor = 20./3.;
                            Delta_factor = 20.;
                        }
                        else{
                            U_factor = 20./3.;
                            Delta_factor = 20.;
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
                            U_factor = 40. / 3.;
                            Delta_factor = 40.;
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
                            U_factor = 20. / 3.;
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

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_tlower_aux() const {return b.t_lower;};
    //double get_tupper_aux() const {return b.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
//
    //void get_freqs_w(double& w, int i) const {w = b.ws[i];};
    //void get_freqs_aux(double& w, int iw) const {w = b.ts[iw];};
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

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
    //auto get_freqGrid_f() const -> FrequencyGrid {return f;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_wlower_f() const {return f.w_lower;};
    //double get_wupper_f() const {return f.w_upper;};
    //double get_tlower_b_aux() const {return b.t_lower;};
    //double get_tupper_b_aux() const {return b.t_upper;};
    //double get_tlower_f_aux() const {return f.t_lower;};
    //double get_tupper_f_aux() const {return f.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_f(double w) const -> double {return f.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    //auto gridtransf_inv_f(double t) const -> double {return f.grid_transf_inv(t);};
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

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
    //auto get_freqGrid_f() const -> FrequencyGrid {return f;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_wlower_f() const {return f.w_lower;};
    //double get_wupper_f() const {return f.w_upper;};
    //double get_tlower_b_aux() const {return b.t_lower;};
    //double get_tupper_b_aux() const {return b.t_upper;};
    //double get_tlower_f_aux() const {return f.t_lower;};
    //double get_tupper_f_aux() const {return f.t_upper;};
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
