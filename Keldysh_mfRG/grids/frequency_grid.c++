#include "frequency_grid.hpp"


/**
 * This function initializes the frequency grid according to the grid parameters
 * w_upper(=-w_lower):  upper limit of frequency box
 * number_of_gridpoints:                 number of frequency points
 * W_scale:             non-linearity of the t_from_frequency function
 */
void FrequencyGrid::initialize_grid() {
    derive_auxiliary_parameters();
    double W;
    for(int i=0; i<number_of_gridpoints; ++i) {
        W = t_lower + i * spacing_auxiliary_gridpoint;
        all_frequencies[i] = frequency_from_t(W);
        assert(isfinite(all_frequencies[i]));
        if (!KELDYSH && !ZERO_T){
            if (type == 'b') all_frequencies[i] = round2bfreq(all_frequencies[i]);
            else             all_frequencies[i] = round2ffreq(all_frequencies[i]);
        }
        auxiliary_grid[i]= t_from_frequency(all_frequencies[i]);
    }
    if (number_of_gridpoints % 2 == 1) {
        all_frequencies[(int) number_of_gridpoints / 2] = 0.;  // make sure that the center of the grid is exactly zero (and not ~10^{-30})
        auxiliary_grid[(int) number_of_gridpoints / 2] = 0.;
    }
    if (!KELDYSH && !ZERO_T) assert (is_doubleOccurencies(all_frequencies) == 0);
}


/// derive auxiliary parameters from essential parameters
void FrequencyGrid::derive_auxiliary_parameters() {
    t_upper = t_from_frequency(w_upper);
    t_lower = t_from_frequency(w_lower);
    spacing_auxiliary_gridpoint = (t_upper - t_lower) / ((double) (number_of_gridpoints-1));
}
void FrequencyGrid::guess_essential_parameters(double Lambda) {
    switch (type) {
        case 'b':
            switch (diag_class) {
                case 1:
                    number_of_gridpoints = nBOS;
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
                    number_of_gridpoints = nBOS2;
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
                    number_of_gridpoints = nBOS3;
                    break;
                default:;
            }
            break;
        case 'f':
            switch (diag_class) {
                case 1:
                    number_of_gridpoints = nFER;
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
                    number_of_gridpoints = nFER2;
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
                    number_of_gridpoints = nFER3;
                    break;
                default:;
            }
            break;
        default:;
    }

    double scale;
    if (REG==2) {
        // scale the grid with Delta until Delta = T_K (Kondo temperature), then scale with T_K
        double Delta = (Lambda + glb_Gamma) / 2.;
        scale = Delta_factor * Delta;
    }
    else if (REG==3) {
        scale = std::max(U_factor * glb_U, Delta_factor * (glb_Gamma) / 2. + Lambda * (Lambda + 1));
    }
    else {
        scale = std::max(U_factor * glb_U, Delta_factor * (glb_Gamma) / 2.);
    }

    all_frequencies = rvec(number_of_gridpoints);
    auxiliary_grid = rvec(number_of_gridpoints);
    set_essential_parameters(scale*15., scale);
    initialize_grid();
}

void FrequencyGrid::update_Wscale(double Wscale) {
    W_scale = Wscale;
    initialize_grid();
}

/** This function returns the index corresponding to the frequency w_in.
 *  It rounds down due to the narrowing conversion from double to int.
 *  This is only used for (linear) interpolations. Hence the narrowing conversion is harmless.
 */
auto FrequencyGrid::get_grid_index(const double w_in) const -> int {
#ifdef PARAMETRIZED_GRID
    double t = t_from_frequency(w_in);

    t = (t - t_lower) / spacing_auxiliary_gridpoint;
#ifdef DENSEGRID
    auto index = ((int) (t + 0.1 )) ;
#else
    auto index = ((int) t) ;
    if (INTERPOLATION==linear) {
        index = std::max(0, index);
        index = std::min(number_of_gridpoints - 2, index);
    }
    else {
        if (all_frequencies[index+1] < w_in) index++;
        index = std::max(0, index);
        index = std::min(number_of_gridpoints - 2, index);
    }
#endif
    return index;

#else
    int j;
    if (INTERPOLATION==linear) {locate(all_frequencies, number_of_gridpoints, w_in, j, 0, number_of_gridpoints-1);} // we cannot interpolate with infinity
    else {locate(all_frequencies, number_of_gridpoints, w_in, j, 0, number_of_gridpoints-1); }
    int index = j;
        assert(all_frequencies[index] - w_in  <= 1e-5 or index == 0);
        assert(w_in - all_frequencies[index+1] < 1e-5 or index == number_of_gridpoints-1);
    return index;
#endif
}

/** This function returns the index corresponding to the frequency w_in.
 *  It rounds down due to the narrowing conversion from double to int.
 *  This is only used for (linear) interpolations. Hence the narrowing conversion is harmless.
 */
auto FrequencyGrid::get_grid_index(double& t, double w_in) const -> int {
    t = t_from_frequency(w_in);
#ifdef PARAMETRIZED_GRID

    double t_rescaled = (t - t_lower) / spacing_auxiliary_gridpoint;
    auto index = ((int) t_rescaled) ;  // round down
    if constexpr(INTERPOLATION==linear) {
        index = std::max(0, index);
        index = std::min(number_of_gridpoints - 2, index);
    }
    else {
        if (all_frequencies[index+1] < w_in and index < number_of_gridpoints-1)
            index++;
        if (all_frequencies[index] > w_in and index > 0)
            index--;
        index = std::max(0, index);
        index = std::min(number_of_gridpoints - 2, index);
        assert(all_frequencies[index] - w_in  <= 1e-5*std::abs(w_in) or index == 0); /// TODO: If this is not satisfied -> use locate
        assert(w_in - all_frequencies[index+1] < 1e-5 or index == number_of_gridpoints-2);
    }
    return index;

#else
    int j;
    if (INTERPOLATION==linear) {locate(all_frequencies, number_of_gridpoints, w_in, j, 0, number_of_gridpoints-1);} // we cannot interpolate with infinity
    else {locate(auxiliary_grid, number_of_gridpoints, t, j, 0, number_of_gridpoints-1); }
    int index = j;
        assert(all_frequencies[index] - w_in  <= 1e-5 or index == 0);
        assert(w_in - all_frequencies[index+1] < 1e-5 or index == number_of_gridpoints-1);
    return index;
#endif

}

/**
 * The used t_from_frequency function is determined here
 * @param w     frequency
 * @return
 */
auto FrequencyGrid::t_from_frequency(double w) const -> double {
    if (KELDYSH) return grid_transf_v2(w, this->W_scale);
    else if (this->type == 'f' and this->diag_class == 1) {
        if (ZERO_T) return grid_transf_v3(w, this->W_scale);
        else return grid_transf_lin(w, this->W_scale);

    }
        //else if (this->type == 'b' and this->diag_class == 1) {
        //    return grid_transf_v2(w, this->W_scale);
        //}
    else {
        if (ZERO_T) return grid_transf_v4(w, this->W_scale);
        else                   return grid_transf_lin(w, this->W_scale);
    }
}

/**
 * The used frequency_from_t function is determined here
 * @param t     point on auxiliary grid
 * @return
 */
auto FrequencyGrid::frequency_from_t(double t) const -> double {
    if (KELDYSH) return grid_transf_inv_v2(t, this->W_scale);
    else if (this->type == 'f' and this->diag_class == 1) {
        if (ZERO_T) return grid_transf_inv_v3(t, this->W_scale);
        else return grid_transf_inv_lin(t, this->W_scale);
    }
        //else if (this->type == 'b' and this->diag_class == 1) {
        //    return grid_transf_inv_v2(w, this->W_scale);
        //}
    else { // TODO(medium): Remove commented part?
        if (KELDYSH || ZERO_T) return grid_transf_inv_v4(t, this->W_scale);
        else return grid_transf_inv_lin(t, this->W_scale);
    }
}

/**
 * This function picks a suitable W_scale for a given wmax for Matsubara T>0
 * @return Wscale    non-linearity of t_from_frequency()
 * @param w1        first positive Matsubara frequency
 * @param wmax      upper bound of frequency grid
 * @param N         relates tmax to t1 by t1=tmax/N (for bosons: (nBOS-1)/2; for fermions: nFER-1)
 */
auto FrequencyGrid::wscale_from_wmax(double & Wscale, const double w1, const double wmax, const int N) -> double {
    if (!KELDYSH && !ZERO_T){
        if (this->type == 'f' and this->diag_class == 1) {
            return wscale_from_wmax_v3(Wscale, w1, wmax, N);
        }
        else {
            return wscale_from_wmax_v1(Wscale, w1, wmax, N);
        }
    }
    else {
        assert(false);
        return 0;
    }
}


/*******************************************    FREQUENCY GRID    *****************************************************/


/**
 * Here are several functions which map the real axis to the compact interval [-1,1]
 * W_scale: sets the scale separating small frequencies (containing interesting structures) and tails
 */

double grid_transf_v1(const double w, const double W_scale) {
    // Version 1: linear around w=0, good for w^(-2) tails
    return w/sqrt(W_scale*W_scale + w*w);
}
double grid_transf_inv_v1(const double t, const double W_scale) {
    // Version 1: linear around w=0, good for w^(-2) tails
    return W_scale*t/sqrt(1.-t*t);
}
double integration_measure_v1(const double t, const double W_scale) {
    double temp = sqrt(1 - t*t);
    return W_scale / (temp*temp*temp);
}


double grid_transf_v2(const double w, const double W_scale) {
    // Version 2: quadratic around w=0, good for w^(-2) tails
    double w2 = w * w;
    return sgn(w) * sqrt((sqrt(w2*w2 + 4 * w2 * W_scale * W_scale) - w2) / 2.) / W_scale;
}
double grid_transf_inv_v2(double t, double W_scale) {
    // Version 2: quadratic around w=0, good for w^(-2) tails
    return W_scale * t * std::abs(t) / sqrt(1. - t * t);
}
double integration_measure_v2(const double t, const double W_scale) {
    double temp = sqrt(1 - t*t);
    return W_scale * sgn(t) * t * (2 - t*t) / (temp*temp*temp);
}
//double grid_transf_v2b(const double w, const double W_scale) {
//    return sgn(w) * log(1 + std::abs(w)/w_a) / log(1. + glb_w_upper/w_a);
//}
//double grid_transf_inv_v2b(const double t, const double W_scale) {
//    return sgn(t) * w_a * (exp(log(1. + glb_w_upper/w_a) * std::abs(t)) - 1);
//}

double grid_transf_v3(const double w, const double W_scale) {
    // Version 3: linear around w=0, good for w^(-1) tails
    const double almost_zero = 1e-12;
    return (std::abs(w) < almost_zero) ? 0. : (-W_scale + sqrt(4*w*w + W_scale*W_scale))/2/w;
}
double grid_transf_inv_v3(const double t, const double W_scale) {
    // Version 3: linear around w=0, good for w^(-1) tails
    return W_scale * t / (1.-t*t);
}
double integration_measure_v3(const double t, const double W_scale) {
    double temp = t*t;
    return W_scale * (1 + temp) / (1 - temp) / (1 - temp);
}

double grid_transf_v4(const double w, const double W_scale) {
    // Version 4: quadratic around w=0, good for w^(-1) tails
    return sgn(w) * sqrt(std::abs(w)/(std::abs(w) + W_scale));
}
double grid_transf_inv_v4(const double t, const double W_scale) {
    // Version 4: quadratic around w=0, good for w^(-1) tails
    return W_scale * sgn(t) * t*t /(1.-t*t);
}
double integration_measure_v4(const double t, const double W_scale) {
    double temp = 1 - t*t;
    return 2 * W_scale * sgn(t) * t / temp / temp;
}



////    linear grid
double grid_transf_lin(double w, double W_scale) {
    return w / W_scale;
}
double grid_transf_inv_lin(double W, double W_scale) {
    return W * W_scale;
}

/**
 * Makes sure that lower bound for W_scale fulfilled
 * Wscale:  current value for W_scale
 * w1:      smallest positive Matsubara frequency (for bosons: 2*M_PI*glb_T, for fermions: M_PI*glb_T)
 * wmax:    maximal frequency
 * N:       relates tmax to t1 by t1=tmax/N (for bosons: (nBOS-1)/2; for fermions: nFER-1)
*/
double wscale_from_wmax_v1(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 1: linear around w=0, good for w^(-2) tails
    Wscale_candidate = wmax * sqrt( (N*N -1*2) / (pow(wmax/w1, 2) - N*N));
    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v2(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 2: quaadratic around w=0, good for w^(-2) tails
    assert(w1 / wmax < 1./(N*N));       // if this fails, then change wmax or use Version 1
    Wscale_candidate = w1 * wmax * N * sqrt(N*N - 1) *sqrt(wmax*wmax - N*N*w1*w1) / std::abs(wmax*wmax - w1*w1*N*N*N*N);

    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v3(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 3: quadratic around w=0, good for w^(-1) tails
    Wscale_candidate = w1 * wmax * (N*N - 1*2.) / sqrt(N * (N * (wmax*wmax + w1*w1) - N*N*w1*wmax - w1*wmax));

    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_lin(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // linear grid
    Wscale_candidate = wmax + 2*M_PI*glb_T;
    return std::max(Wscale, Wscale_candidate);
}






/*******************************************    HYBRID GRID    *****************************************************/

void hybridGrid::derive_auxiliary_parameters() {
    const double temp = (pos_section_boundaries[0] + 2.*pos_section_boundaries[1]);
    const double temp_quad = temp * w_upper - pos_section_boundaries[1]*pos_section_boundaries[1];
    aux_pos_section_boundaries[0] = (2.*t_upper*pos_section_boundaries[0]*w_upper)                              / temp_quad;
    aux_pos_section_boundaries[1] = (t_upper*(pos_section_boundaries[0] + pos_section_boundaries[1]) * w_upper) / temp_quad;

    recip_curvature_quad = (4.*t_upper*t_upper*pos_section_boundaries[0]*w_upper*w_upper) / (temp_quad*temp_quad );
    recip_slope_lin = t_upper*w_upper / (temp_quad);
    factor_rat = pos_section_boundaries[1]*pos_section_boundaries[1] / (temp);
    rescale_rat = (t_upper * w_upper * temp) / temp_quad;
    spacing_auxiliary_gridpoint = 2.*t_upper / ((double)number_of_gridpoints - 1.);
}

void hybridGrid::guess_essential_parameters(const double Lambda) {
    const double Delta = (glb_Gamma + Lambda) * 0.5;

    switch (diag_class) {
        case 1:
            if (type == 'b') {
                number_of_gridpoints = nBOS;
                pos_section_boundaries[0] = 20.*Delta;
                pos_section_boundaries[1] = 30.*Delta;
                w_upper = 40.*15.*Delta;
            }
            else {
                assert(type == 'f');
                number_of_gridpoints = nFER;
                pos_section_boundaries[0] = 10.;
                pos_section_boundaries[1] = 50.*Delta;
                w_upper = 40.*15.*Delta;
            }
            break;
        case 2:
            if (type == 'b') {
                number_of_gridpoints = nBOS2;
                pos_section_boundaries[0] = Delta;
                pos_section_boundaries[1] = 3.*Delta;
                w_upper = 40.*15.*Delta;
            }
            else {
                assert(type == 'f');
                number_of_gridpoints = nFER2;
                pos_section_boundaries[0] = Delta;
                pos_section_boundaries[1] = 3.*Delta;
                w_upper = 20.*15.*Delta;
            }
            break;
        case 3:
            if (type == 'b') {
                number_of_gridpoints = nBOS3;
                pos_section_boundaries[0] = Delta;
                pos_section_boundaries[1] = 3.*Delta;
                w_upper = 30.*15.*Delta;
            }
            else {
                assert(type == 'f');
                number_of_gridpoints = nFER3;
                pos_section_boundaries[0] = Delta;
                pos_section_boundaries[1] = 3.*Delta;
                w_upper = 15*15.*Delta;
            }
            break;
        default:
            break;
    }

    w_lower = -w_upper;
    all_frequencies = rvec(number_of_gridpoints);
    auxiliary_grid = rvec(number_of_gridpoints);
    initialize_grid();
}

double hybridGrid::frequency_from_t(const double t) const {
    const double t_abs = std::abs(t);
    if (t_abs < aux_pos_section_boundaries[0]) {
        // quaddratic part
        const double result = t_abs*t_abs / recip_curvature_quad * sgn(t);
        assert(isfinite(result));
        return result;
    }
    else if (t_abs < aux_pos_section_boundaries[1]) {
        // linear part
        const double result = (pos_section_boundaries[0] + (t_abs - aux_pos_section_boundaries[0]) / recip_slope_lin)* sgn(t);
        assert(isfinite(result));
        return result;
    }
    else {
        // rational part
        const double result = factor_rat / (1. - t_abs / rescale_rat) * sgn(t);
        assert(isfinite(result));
        return result;
    }

}

double hybridGrid::t_from_frequency(const double w) const {
    const double w_abs = std::abs(w);
    if (w_abs < pos_section_boundaries[0]) {
        // quadratic part
        const double result = sqrt(w_abs * recip_curvature_quad) * sgn(w);
        assert(isfinite(result));
        return result;
    }
    else if (w_abs < pos_section_boundaries[1]) {
        // linear part
        const double result = ((w_abs - pos_section_boundaries[0]) * recip_slope_lin + aux_pos_section_boundaries[0]) * sgn(w);
        assert(isfinite(result));
        return result;
    }
    else {
        // rational part
        const double result = (1. - factor_rat / w_abs) * rescale_rat * sgn(w);
        assert(isfinite(result));
        return result;
    }
}

void hybridGrid::initialize_grid() {
    derive_auxiliary_parameters();
    for (int i = 0; i < number_of_gridpoints; i++) {
        const double t = t_lower + spacing_auxiliary_gridpoint * (double)i;
        auxiliary_grid[i] = t;
        all_frequencies[i] = frequency_from_t(t);
        if constexpr(!KELDYSH && !ZERO_T){
            if (type == 'b') all_frequencies[i] = round2bfreq(all_frequencies[i]);
            else             all_frequencies[i] = round2ffreq(all_frequencies[i]);
        }
    }
    if (number_of_gridpoints % 2 == 1) {
        all_frequencies[(int) number_of_gridpoints / 2] = 0.;  // make sure that the center of the grid is exactly zero (and not ~10^{-30})
        auxiliary_grid[(int) number_of_gridpoints / 2] = 0.;
    }
    if (!KELDYSH && !ZERO_T) assert (is_doubleOccurencies(all_frequencies) == 0);
}

int hybridGrid::get_grid_index(const double frequency)const {
    /// Do I need this? The next function basically does the same
    const double t = t_from_frequency(frequency);
#ifdef DENSEGRID
    auto index = ((int) ((t - t_base) / spacing_auxiliary_gridpoint + 0.1 )) ;
#else
    int index = int((t - t_lower) / spacing_auxiliary_gridpoint + 1e-12);
    index = std::max(0, index);
    index = std::min(number_of_gridpoints - 2, index);
#endif
    assert(all_frequencies[index] <= frequency);
    assert(all_frequencies[index+1] >= frequency);
    assert(auxiliary_grid[index] <= t);
    assert(auxiliary_grid[index+1] >= t);
    return index;
}

int hybridGrid::get_grid_index(double& t, const double frequency) const{
    t = t_from_frequency(frequency);
#ifdef DENSEGRID
    auto index = ((int) ((t - t_base) / spacing_auxiliary_gridpoint + 0.1 )) ;
#else
    int index = int((t - t_lower) / spacing_auxiliary_gridpoint + 1e-12);
    index = std::max(0, index);
    index = std::min(number_of_gridpoints - 2, index);
#endif
    assert(all_frequencies[index] <= frequency);
    assert(all_frequencies[index+1] >= frequency);
    assert(auxiliary_grid[index] <= t);
    assert(auxiliary_grid[index+1] >= t);
    return index;
}