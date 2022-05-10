#ifndef KELDYSH_MFRG_BUBBLE_HPP
#define KELDYSH_MFRG_BUBBLE_HPP

#include <cmath>                            // for using the macro M_PI as pi
#include "../symmetries/Keldysh_symmetries.hpp"  // for independent Keldysh components and utilities
#include "../correlation_functions/four_point/vertex.hpp"                         // vertex class
#include "../correlation_functions/two_point/selfenergy.hpp"                     // self-energy class
#include "../correlation_functions/two_point/propagator.hpp"                     // propagator class
#include "../integrator/integrator.hpp"          // integration routines
#include "../utilities/util.hpp"                 // measuring time, printing text output
#include "../utilities/mpi_setup.hpp"            // mpi parallelization routines
#include "../asymptotic_corrections/correction_functions.hpp"            // correction terms due to finite integration range
#include "../utilities/write_data2file.hpp"      // write vectors into hdf5 file
#include "../grids/momentum_grid.hpp"            // Momentum grid specific to the 2D Hubbard model


/// Class combining two propagators, either GG or GS+SG
template <typename Q>
class Bubble{
public:
    const Propagator<Q>& g; // Access needed when computing the full bubble
    const Propagator<Q>& s;
    const bool diff;
    /**
     * Constructor:
     * @param propagatorG : first propagator (always a standard one)
     * @param propagatorS : second propagator (standard or single-scale/differentiated, depending on "diff_in")
     * @param diff_in      : whether to compute standard (false) or differentiated (true) bubble
     */
    Bubble(const Propagator<Q>& propagatorG, const Propagator<Q>& propagatorS, const bool diff_in)
        :g(propagatorG), s(propagatorS), diff(diff_in) {
        g.initInterpolator();
        s.initInterpolator();
    };

    /**
     * Call operator:
     * @param iK    : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param v1    : frequency of first propagator
     * @param v2    : frequency of second propagator
     * @param i_in  : internal structure index
     * @return comp : value of the bubble evaluated at (iK, v1, v2)
     */
    auto value(int iK, double v1, double v2, int i_in) const -> Q{
        Q ans;
        if(diff){
            if (KELDYSH){
                switch (iK) {
                    case 3: //AA
                        ans = myconj(g.valsmooth(0, v1, i_in)) * myconj(s.valsmooth(0, v2, i_in)) + myconj(s.valsmooth(0, v1, i_in)) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 6: //AR
                        ans = myconj(g.valsmooth(0, v1, i_in)) * s.valsmooth(0, v2, i_in) + myconj(s.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                        break;
                    case 7: //AK
                        ans = myconj(g.valsmooth(0, v1, i_in)) * s.valsmooth(1, v2, i_in) + myconj(s.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                        break;
                    case 9: //RA
                        ans = g.valsmooth(0, v1, i_in) * myconj(s.valsmooth(0, v2, i_in)) + s.valsmooth(0, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 11://KA
                        ans = g.valsmooth(1, v1, i_in) * myconj(s.valsmooth(0, v2, i_in)) + s.valsmooth(1, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 12://RR
                        ans = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                        break;
                    case 13://RK
                        ans = g.valsmooth(0, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    case 14://KR
                        ans = g.valsmooth(1, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                        break;
                    case 15://KK
                        ans = g.valsmooth(1, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(1, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    default:
                        return 0.;
                }
            }
            else{ // Matsubara
                ans = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
            }
        }
        else {
            if (KELDYSH){
                switch (iK){ // labelling propagators from top (t: left) to bottom (t: right); a,t: G(v+w/2)G(v-w/2), p: G(w/2-v)G(w/2+v)
                    case 3: //AA
                        ans = myconj(g.valsmooth(0, v1, i_in)) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 6: //AR
                        ans = myconj(g.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                        break;
                    case 7: //AK
                        ans = myconj(g.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                        break;
                    case 9: //RA
                        ans = g.valsmooth(0, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 11://KA
                        ans = g.valsmooth(1, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 12://RR
                        ans = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                        break;
                    case 13://RK
                        ans = g.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    case 14://KR
                        ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                        break;
                    case 15://KK
                        ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(1, v2, i_in);
                        break;
                    default:
                        return 0.;
                }
            }
            else{ // Matsubara
                ans = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
            }
        }
        if (!KELDYSH) {
            if (PARTICLE_HOLE_SYMMETRY) {
                ans *= -1.;     // -1=glb_i^2; needed for particle-hole symmetry in Matsubara (we only save the imaginary part of self-energy and propagators)
            }
        }

        assert(isfinite(ans) == true);
        return ans;
    }

    template<char ch_bubble>
    Q value_Matsubara(const double w, const double vpp, int i_in) const {
        Q result;

        double v1 = convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w);
        double v2 = convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w);
        if (diff) {
            result = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);

        }
        else {
            result = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);

        }
        if constexpr(PARTICLE_HOLE_SYMMETRY) {
            result *= -1.;     // -1=glb_i^2; needed for particle-hole symmetry in Matsubara (we only save the imaginary part of self-energy and propagators)
        }
        assert(isfinite(result));
        return result;
    }

    /**
     * Wrapper for value function above, providing the natural arguments for evaluation of the bubble in each channel:
     * @param iK      : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param w       : bubble transfer frequency of the corresponding channel
     * @param vpp     : bubble integration frequency of the corresponding channel
     * @param i_in    : internal structure index
     * @param channel : channel to which the bubble belongs
     * @return Q      : value of the bubble evaluated at the arguments described above (usually comp)
     */
    auto value(int iK, double w, double vpp, int i_in, char channel) const -> Q {
        double wa_1, wa_2, wp_1, wp_2, wt_1, wt_2;
        if (KELDYSH || ZERO_T){
            wa_1 = wa_2 = wp_1 = wp_2 = wt_1 = wt_2 = w / 2.;
        }
        else{ // bosonic frequencies have to be rounded to an integer value for finite-temperature Matsubara calculation
            wa_1 = wp_2 = wt_1 = floor2bfreq(w / 2.);
            wa_2 = wp_1 = wt_2 = ceil2bfreq(w / 2.);
        } // TODO(medium): Put this first part into an extra function?
        Q Pival;
        switch (channel) {
            case 'a':
                Pival = value(iK, vpp - wa_1, vpp + wa_2, i_in);    //vppa-1/2wa, vppa+1/2wa for the a-channel
                break;
            case 'p':
                Pival = value(iK, wp_1 + vpp, wp_2 - vpp, i_in);    //wp/2+vppp, wp/2-vppp for the p-channel
                break;
            case 't':
                Pival = value(iK, vpp - wt_1, vpp + wt_2, i_in);    //vppt-1/2wt, vppt+1/2wt for the t-channel
                break;
            default:;
        }
        assert(isfinite(Pival) == true);
        return Pival;
    }

    template<char ch_bubble>
    auto convert_to_fermionic_frequencies_1(double vpp, double w) const -> double {
        if constexpr(KELDYSH || ZERO_T) {
            if constexpr(ch_bubble == 'p') return vpp + w * 0.5;
            else return vpp - w * 0.5;
        }
        else { // Matsubata zeroT
            if constexpr(ch_bubble == 'p') return vpp + ceil2bfreq(w * 0.5);
            else return vpp - floor2bfreq(w * 0.5);
        }
    }
    template<char ch_bubble>
    auto convert_to_fermionic_frequencies_2(double vpp, double w) const -> double {
        if constexpr(KELDYSH || ZERO_T) {
            if constexpr(ch_bubble == 'p') return w * 0.5 - vpp;
            else return vpp + w * 0.5;
        }
        else { // Matsubata zeroT
            if constexpr(ch_bubble == 'p') return floor2bfreq(w * 0.5) - vpp;
            else return vpp + ceil2bfreq(w * 0.5);
        }
    }


    template<char ch_bubble>
    auto value_vectorized(const double w, const double vpp, const int i_in) const {
        if constexpr(KELDYSH) {
            //static_assert(KELDYSH, "vector-valued integrand only allowed for Keldysh formalism.");
            using result_type = Eigen::Matrix<Q, 4, 4>;
            using buffertype_propagator = Eigen::Matrix<Q, 2, 2>;
            const buffertype_propagator G1 = g.template valsmooth_vectorized<buffertype_propagator>(convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
            const buffertype_propagator G2 = g.template valsmooth_vectorized<buffertype_propagator>(convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);

            const Q G00_1 = G1(0,0); // g.valsmooth(0, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
            const Q G01_1 = G1(0,1); // g.valsmooth(1, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
            const Q G10_1 = G1(1,0); // g.valsmooth(2, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
            const Q G11_1 = G1(1,1); // g.valsmooth(3, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
            const Q G00_2 = G2(0,0); // g.valsmooth(0, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
            const Q G01_2 = G2(0,1); // g.valsmooth(1, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
            const Q G10_2 = G2(1,0); // g.valsmooth(2, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
            const Q G11_2 = G2(1,1); // g.valsmooth(3, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);

            result_type result;

            if (diff) {
                const buffertype_propagator S1 = s.template valsmooth_vectorized<buffertype_propagator>(convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);;
                const buffertype_propagator S2 = s.template valsmooth_vectorized<buffertype_propagator>(convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);;


                const Q S00_1 = S1(0,0); // g.valsmooth(0, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
                const Q S01_1 = S1(0,1); // g.valsmooth(1, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
                const Q S10_1 = S1(1,0); // g.valsmooth(2, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
                const Q S11_1 = S1(1,1); // g.valsmooth(3, convert_to_fermionic_frequencies_1<ch_bubble>(vpp, w), i_in);
                const Q S00_2 = S2(0,0); // g.valsmooth(0, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
                const Q S01_2 = S2(0,1); // g.valsmooth(1, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
                const Q S10_2 = S2(1,0); // g.valsmooth(2, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);
                const Q S11_2 = S2(1,1); // g.valsmooth(3, convert_to_fermionic_frequencies_2<ch_bubble>(vpp, w), i_in);

                if constexpr(ch_bubble=='p') {
                    result <<    S00_1*G00_2 + G00_1*S00_2, S00_1*G01_2 + G00_1*S01_2, S01_1*G00_2 + G01_1*S00_2, S01_1*G01_2 + G01_1*S01_2,
                                 S00_1*G10_2 + G00_1*S10_2, S00_1*G11_2 + G00_1*S11_2, S01_1*G10_2 + G01_1*S10_2, S01_1*G11_2 + G01_1*S11_2,
                                 S10_1*G00_2 + G10_1*S00_2, S10_1*G01_2 + G10_1*S01_2, S11_1*G00_2 + G11_1*S00_2, S11_1*G01_2 + G11_1*S01_2,
                                 S10_1*G10_2 + G10_1*S10_2, S10_1*G11_2 + G10_1*S11_2, S11_1*G10_2 + G11_1*S10_2, S11_1*G11_2 + G11_1*S11_2;
                }
                else { // a or t -channel
                    result <<    S00_1*G00_2 + G00_1*S00_2, S00_1*G10_2 + G00_1*S10_2, S01_1*G00_2 + G01_1*S00_2, S01_1*G10_2 + G01_1*S10_2,
                                 S00_1*G01_2 + G00_1*S01_2, S00_1*G11_2 + G00_1*S11_2, S01_1*G01_2 + G01_1*S01_2, S01_1*G11_2 + G01_1*S11_2,
                                 S10_1*G00_2 + G10_1*S00_2, S10_1*G10_2 + G10_1*S10_2, S11_1*G00_2 + G11_1*S00_2, S11_1*G10_2 + G11_1*S10_2,
                                 S10_1*G01_2 + G10_1*S01_2, S10_1*G11_2 + G10_1*S11_2, S11_1*G01_2 + G11_1*S01_2, S11_1*G11_2 + G11_1*S11_2;
                }
            }
            else {
                if constexpr(ch_bubble=='p') {
                    result <<    G00_1*G00_2, G00_1*G01_2, G01_1*G00_2, G01_1*G01_2,
                                 G00_1*G10_2, G00_1*G11_2, G01_1*G10_2, G01_1*G11_2,
                                 G10_1*G00_2, G10_1*G01_2, G11_1*G00_2, G11_1*G01_2,
                                 G10_1*G10_2, G10_1*G11_2, G11_1*G10_2, G11_1*G11_2;
                }
                else { // a or t -channel
                    result <<    G00_1*G00_2, G00_1*G10_2, G01_1*G00_2, G01_1*G10_2,
                                 G00_1*G01_2, G00_1*G11_2, G01_1*G01_2, G01_1*G11_2,
                                 G10_1*G00_2, G10_1*G10_2, G11_1*G00_2, G11_1*G10_2,
                                 G10_1*G01_2, G10_1*G11_2, G11_1*G01_2, G11_1*G11_2;
                }
            }
            return result;

        }
        else {
            Q result = value_Matsubara<ch_bubble>(w, vpp, i_in);
            return result;
        }

    }
};

#endif //KELDYSH_MFRG_BUBBLE_HPP
