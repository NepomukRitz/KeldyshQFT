#ifndef KELDYSH_MFRG_INTEGRANDSE_HPP
#define KELDYSH_MFRG_INTEGRANDSE_HPP

#include <cmath>                    // for using the macro M_PI as pi
#include "../correlation_functions/two_point/selfenergy.hpp"             // self-energy class
#include "../correlation_functions/four_point/vertex.hpp"                 // vertex class
#include "../correlation_functions/two_point/propagator.hpp"             // propagator class
#include "../parameters/master_parameters.hpp"             // system parameters (vector lengths etc.)
#include "../integrator/integrator.hpp"             // integration routines
#include "../utilities/write_data2file.hpp"        // save integrand for debugging purposes
#include "../asymptotic_corrections/correction_functions.hpp"    // analytical results for the tails of the loop integral
#include "../utilities/hdf5_routines.hpp"

/// possible tests: ---> See bubble integrand


namespace selfenergy_loop {
    template<bool version> constexpr char pick_channel() {return version == 0 ? 't' : 'a';}
    template<bool version, int ispin> constexpr int pick_spin() {return version == 0 ? ispin : 2-ispin;}
}

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 */

/**
 * Class for the integrand of a loop calculation. Invoked by LoopCalculator.
 * @tparam Q Type of the data.
 * @tparam vertType Type of the vertex to be used.
 * @tparam all_spins Compute all spin components?
 * @tparam return_type Type of the data of the result. Typically = Q.
 * @tparam version Version of the self-energy loop. See documentation for LoopFunctionCalculator.
 */
template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
class IntegrandSE {
#if KELDYSH_FORMALISM
    using buffertype_propagator = Eigen::Matrix<Q,1,4>;
    using buffertype_vertex = Eigen::Matrix<Q,4,myColsAtCompileTime<return_type>()>;
#else
    using buffertype_propagator = Q;
    using buffertype_vertex = Q;
#endif

    const int type;
    std::vector<int> components = std::vector<int>((CONTOUR_BASIS != 1 ? 6 : 16));

    const vertType& vertex;
    const Propagator<Q>& propagator;

    const int iK;
    const double v;
    const int i_in;
    const int i_spin;
    const int i0_vertex_left;
    const int i0_vertex_right;
    Eigen::Matrix<Q,4, 4> transform_to_KeldyshBasis;

    void set_Keldysh_components_to_be_calculated();

    return_type Keldysh_value(double vp) const;
    Q Matsubara_value(double vp) const;

    void evaluate_propagator(Q& Gi, const int iK, const double vp) const;
    void evaluate_propagator(Q& GM, const double vp) const; // Matsubara version
    auto evaluate_propagator_vectorized(const double vp) const -> buffertype_propagator;

    void evaluate_vertex(Q& factorClosedAbove, Q& factorAClosedBelow, const int iK, const double vp) const; // for symmetrized Keldysh flow
    void evaluate_vertex(Q& factorClosedAbove, const int iK, const double vp) const; // for unsymmetrized Keldysh/Matsubara flow
    auto evaluate_vertex_vectorized(const double vp) const -> buffertype_vertex;

public:

    /**
     *
     * @param type_in determines Keldysh components (always type = 0 for MF) is the element in Sigma matrix to be computed:
     * ```
     *    0  1    ___     K  A
     *    2  3            R  0
     * ```
     * For example: \n
     *      type = 0 -> Keldysh component \n
     *      type = 2 -> Retarded component
     * @param vertex_in Reference to the input vertex.
     * @param prop_in Reference to the input propagator.
     * @param iK_in Keldysh index
     * @param i_spin_in spin index
     * @param v_in external fermionic frequency index
     * @param i_in_in internal index (trivially 0 most of the time)
     */
    IntegrandSE(const int type_in, const vertType& vertex_in, const Propagator<Q>& prop_in,
                const int iK_in, const int i_spin_in, const double v_in, const int i_in_in)
                :type(type_in), vertex(vertex_in), propagator(prop_in), iK(iK_in), v(v_in), i_in(i_in_in), i_spin(i_spin_in),
                 i0_vertex_left(type_in*4), i0_vertex_right(type_in*4){
        if constexpr(KELDYSH){
            set_Keldysh_components_to_be_calculated();
            //if constexpr(DEBUG_SYMMETRIES) {
            transform_to_KeldyshBasis << 0.5, 0.5, 0.5, 0.5,
                                        -0.5, 0.5,-0.5, 0.5,
                                        -0.5,-0.5, 0.5, 0.5,
                                         0.5,-0.5,-0.5, 0.5;

        }
    }

    /**
     * Call operator:
     * @param vp frequency at which to evaluate integrand (to be integrated over)
     * @return value of the integrand object evaluated at frequency vp
     */
    auto operator()(double vp) const -> return_type;

    void save_integrand() const;
    void save_integrand(const rvec& freqs, const std::string& filename_prefix) const;
    void get_integrand_vals(const rvec& freqs, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& integrand_vals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals)  const;

};

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::set_Keldysh_components_to_be_calculated() {
#if not SWITCH_SUM_N_INTEGRAL
    if constexpr(version == 0) {
        #if CONTOUR_BASIS != 1
            if(type==2){  //Check which kind of contribution is calculated
                // retarded selfenergy
                components[0]=3;    //Vertex component associated to Retarded propagator
                components[1]=6;    //Vertex component associated to Advanced propagator
                components[2]=7;    //Vertex component associated to Keldysh propagator
                components[3]=3;    //Vertex component associated to Retarded propagator in symmetrized flow
                components[4]=9;    //Vertex component associated to Advanced propagator in symmetrized flow
                components[5]=11;   //Vertex component associated to Keldysh propagator in symmetrized flow
            }
            else {
                assert(type==0);
                // Keldysh selfenergy
                components[0]=1;    //Vertex component associated to Retarded propagator
                components[1]=4;    //Vertex component associated to Advanced propagator
                components[2]=5;    //Vertex component associated to Keldysh propagator
                components[3]=2;    //Vertex component associated to Retarded propagator in symmetrized flow
                components[4]=8;    //Vertex component associated to Advanced propagator in symmetrized flow
                components[5]=10;   //Vertex component associated to Keldysh propagator in symmetrized flow
            }
        #else
            if(type==0){  //Check which kind of contribution is calculated
                components[0]=0;    //Vertex component associated to propagator with contour indices 00
                components[1]=4;    //Vertex component associated to propagator with contour indices 01
                components[2]=1;    //Vertex component associated to propagator with contour indices 10
                components[3]=5;    //Vertex component associated to propagator with contour indices 11
            }
            else if (type==1){
                components[0]=2;    //Vertex component associated to propagator with contour indices 00
                components[1]=6;    //Vertex component associated to propagator with contour indices 01
                components[2]=3;    //Vertex component associated to propagator with contour indices 10
                components[3]=7;    //Vertex component associated to propagator with contour indices 11
            }
            else if (type==2){
                components[0]=8;    //Vertex component associated to propagator with contour indices 00
                components[1]=12;    //Vertex component associated to propagator with contour indices 01
                components[2]=9;   //Vertex component associated to propagator with contour indices 10
                components[3]=13;   //Vertex component associated to propagator with contour indices 11
            }
            else {
                assert(type == 3);
                components[0]=10;   //Vertex component associated to propagator with contour indices 00
                components[1]=14;   //Vertex component associated to propagator with contour indices 01
                components[2]=11;   //Vertex component associated to propagator with contour indices 10
                components[3]=15;   //Vertex component associated to propagator with contour indices 11
            }
        #endif
    }
    else { // version == 1
        #if CONTOUR_BASIS != 1
            if(type==2){  //Check which kind of contribution is calculated
                // retarded selfenergy
                components[0]=3;    //Vertex component associated to Retarded propagator
                components[1]=5;    //Vertex component associated to Advanced propagator
                components[2]=7;    //Vertex component associated to Keldysh propagator
                components[3]=3;    //Vertex component associated to Retarded propagator in symmetrized flow
                components[4]=10;    //Vertex component associated to Advanced propagator in symmetrized flow
                components[5]=11;   //Vertex component associated to Keldysh propagator in symmetrized flow
            }
            else {
                assert(type==0);
                // Keldysh selfenergy
                components[0]=2;    //Vertex component associated to Retarded propagator
                components[1]=4;    //Vertex component associated to Advanced propagator
                components[2]=6;    //Vertex component associated to Keldysh propagator
                components[3]=1;    //Vertex component associated to Retarded propagator in symmetrized flow
                components[4]=8;    //Vertex component associated to Advanced propagator in symmetrized flow
                components[5]=9;   //Vertex component associated to Keldysh propagator in symmetrized flow
            }
        #else
            if(type==0){  //Check which kind of contribution is calculated
                components[0]=0;    //Vertex component associated to propagator with contour indices 00
                components[1]=4;    //Vertex component associated to propagator with contour indices 01
                components[2]=2;    //Vertex component associated to propagator with contour indices 10
                components[3]=6;    //Vertex component associated to propagator with contour indices 11
            }
            else if (type==1){
                components[0]=1;    //Vertex component associated to propagator with contour indices 00
                components[1]=5;    //Vertex component associated to propagator with contour indices 01
                components[2]=3;    //Vertex component associated to propagator with contour indices 10
                components[3]=7;    //Vertex component associated to propagator with contour indices 11
            }
            else if (type==2){
                components[0]=8;    //Vertex component associated to propagator with contour indices 00
                components[1]=12;    //Vertex component associated to propagator with contour indices 01
                components[2]=10;   //Vertex component associated to propagator with contour indices 10
                components[3]=14;   //Vertex component associated to propagator with contour indices 11
            }
            else {
                assert(type == 3);
                components[0]=9;   //Vertex component associated to propagator with contour indices 00
                components[1]=13;   //Vertex component associated to propagator with contour indices 01
                components[2]=11;   //Vertex component associated to propagator with contour indices 10
                components[3]=15;   //Vertex component associated to propagator with contour indices 11
            }
        #endif
    }
#endif
}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
auto IntegrandSE<Q,vertType,all_spins,return_type,version>::operator()(const double vp) const -> return_type {
    if constexpr(KELDYSH){return Keldysh_value(vp);}
    else{return Matsubara_value(vp);}
}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::get_integrand_vals(const rvec& freqs, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& integrand_vals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals) const {
    int npoints = freqs.size();

    for (int i=0; i<npoints; ++i) {

        double vpp = freqs[i];


        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> integrand_value(1, 4);
        //Q integrand_value;

        integrand_value = (*this)(vpp).transpose();
        integrand_vals.col(i) = integrand_value;

        buffertype_vertex vertex_val = evaluate_vertex_vectorized(vpp);

        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_val_temp = vertex_val;
        vertex_val_temp.resize(16,1);
        vertex_vals.col(i) = vertex_val_temp;
    }


}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = 1e5;

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;

        wl = propagator.selfenergy.Sigma.frequencies.primary_grid.w_lower * 2.;
        wu = propagator.selfenergy.Sigma.frequencies.primary_grid.w_upper * 2.;

        double vpp = wl + i * (wu - wl) / (npoints - 1);
        freqs[i] = vpp;
    }

    save_integrand(freqs, "");

}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::save_integrand(const rvec& freqs, const std::string& filename_prefix) const {
    int npoints = freqs.size();

    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> integrand_vals(4, npoints);
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_vals(VECTORIZED_INTEGRATION and KELDYSH_FORMALISM ? 16 : 1, npoints);

    get_integrand_vals(freqs, integrand_vals, vertex_vals);

    std::string filename = data_dir + filename_prefix + "integrand_SE";
    filename += //"_i0=" + std::to_string(i0)       /// TODO: add this when Elias interchanged order of integration and Keldysh sum
                //+ "_i2=" + std::to_string(i2)
                + "_v=" + std::to_string(v);
    filename += + ".h5";

    utils::print("saving integrand to file ", filename, "\n");
    if (mpi_world_rank() == 0) {
        H5::H5File file = H5::H5File(filename, H5F_ACC_TRUNC);
        write_to_hdf(file, "v", freqs, false);
        write_to_hdf(file, "integrand", integrand_vals, false);
        write_to_hdf(file, "vertex", vertex_vals, false);
    }
}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
auto IntegrandSE<Q,vertType,all_spins,return_type,version>::Keldysh_value(const double vp) const -> return_type {
#if SWITCH_SUM_N_INTEGRAL

    buffertype_propagator G1 = evaluate_propagator_vectorized( vp);
    buffertype_propagator G2 = evaluate_propagator_vectorized(-vp);

    buffertype_vertex V1 = evaluate_vertex_vectorized( vp);
    buffertype_vertex V2 = evaluate_vertex_vectorized(-vp);

if constexpr(std::is_same_v<return_type,double> or std::is_same_v<return_type,comp>) {
    const return_type result = (G1*V1 + G2*V2).eval()[0] * 0.5;
    return result;
}
else {
    if constexpr(CONTOUR_BASIS) {
        const return_type result = (G1*V1*transform_to_KeldyshBasis + G2*V2*transform_to_KeldyshBasis).eval() * 0.5;
        return result;

    }
    else {
        const return_type result = (G1*V1 + G2*V2).eval() * 0.5;
        return result;
    }
}

#else
    Q Gi, Gi2;
    evaluate_propagator(Gi , iK, vp);
    evaluate_propagator(Gi2, iK,-vp);

    Q factorClosedAbove, factorClosedAbove2;
#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
    Q factorClosedBelow;
    evaluate_vertex(factorClosedAbove, factorClosedBelow, iK, vp);
    return (1./2.) * Gi * (factorClosedAbove + factorClosedBelow);
#else
    evaluate_vertex(factorClosedAbove,  iK, vp);
    evaluate_vertex(factorClosedAbove2, iK,-vp);
    return (Gi * factorClosedAbove +  Gi2 * factorClosedAbove2) * 0.5;
#endif  // SYMMETRIZED_SELF_ENERGY_FLOW
#endif  // SWITCH_SUM_N_INTEGRAL

}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
auto IntegrandSE<Q,vertType,all_spins,return_type,version>::Matsubara_value(const double vp) const -> Q {
    Q GM;
    evaluate_propagator(GM, vp);

    Q factorClosedAbove = evaluate_vertex_vectorized(vp);

    if (!PARTICLE_HOLE_SYMMETRY) {
        return (GM * factorClosedAbove);
    }
    else {
        // in the particle-hole symmetric_full case in Matsubara we only save the imaginary part of the selfenergy Im(Sigma)
        // Accordingly the saved propagator is -Im(G)
        // Hence we need an additional factor of -1
        return -(GM * factorClosedAbove);
    }
}

template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_propagator(Q &Gi, const int iK, const double vp) const {
#if CONTOUR_BASIS != 1
    switch (iK) {
        case 0:
            Gi = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
            break;
        case 1:
            Gi = myconj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
            break;
        case 2:
            Gi = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)
            break;
        default:;
    }
#else
    Gi = propagator.valsmooth(iK, vp, i_in);
#endif

}

template<typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_propagator(Q &GM, const double vp) const {
    GM = propagator.valsmooth(0, vp, i_in);           // Matsubara propagator (full or single scale)
}

template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
auto IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_propagator_vectorized(const double vp) const -> buffertype_propagator {
    buffertype_propagator G = propagator.template valsmooth_vectorized<buffertype_propagator>(vp, i_in);
    std::swap(G(1), G(2));
    return G;
}


template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_vertex(Q &factorClosedAbove, Q &factorClosedBelow,
                                     const int iK, const double vp) const {
    using namespace selfenergy_loop;

    VertexInput inputClosedAbove (components[iK]  , i_spin, 0., vp, v, i_in, pick_channel<version>());
    VertexInput inputClosedBelow (components[iK+3], i_spin, 0., v, vp, i_in, pick_channel<version>());
    factorClosedAbove = vertex.template value<pick_channel<version>()>(inputClosedAbove);
    factorClosedBelow = vertex.template value<pick_channel<version>()>(inputClosedBelow);
}

template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
void IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_vertex(Q &factorClosedAbove, const int iK, const double vp) const {
    using namespace selfenergy_loop;
    // "components" are all zero in Matsubara case -> this function also works for Matsubara
    VertexInput inputClosedAbove (components[iK], i_spin, 0, vp, v, i_in, pick_channel<version>());
    factorClosedAbove = vertex.template value<pick_channel<version>()>(inputClosedAbove);
}

template <typename Q, typename vertType, bool all_spins, typename return_type, bool version>
auto IntegrandSE<Q,vertType,all_spins,return_type,version>::evaluate_vertex_vectorized(const double vp) const -> buffertype_vertex {
    using namespace selfenergy_loop;

    if constexpr(all_spins) {
        const VertexInput inputClosedAbove_spin0(i0_vertex_right, 0, 0, vp, v, i_in, pick_channel<version>());
        const VertexInput inputClosedAbove_spin1(i0_vertex_right, 0, 0, vp, v, i_in, pick_channel<version>());
#if VECTORIZED_INTEGRATION == 1 and KELDYSH_FORMALISM
        using valuetype_fetch = Eigen::Matrix<Q,4*return_type::ColsAtCompileTime,1>;
        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> result0 = vertex.template value_symmetry_expanded<pick_spin<version,0>() ,pick_channel<version>(),valuetype_fetch>(inputClosedAbove_spin0);
        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> result1 = vertex.template value_symmetry_expanded<pick_spin<version,1>(),pick_channel<version>(),valuetype_fetch>(inputClosedAbove_spin1);
        result0.resize(4,return_type::ColsAtCompileTime);
        result1.resize(4,return_type::ColsAtCompileTime);
        //static_assert(decltype(result0)::RowsAtCompileTime == 4);
#else
        const buffertype_vertex result0 = vertex.template value_symmetry_expanded<pick_spin<version,0>(),pick_channel<version>(),buffertype_vertex>(inputClosedAbove_spin0);
        const buffertype_vertex result1 = vertex.template value_symmetry_expanded<pick_spin<version,1>(),pick_channel<version>(),buffertype_vertex>(inputClosedAbove_spin1);

#endif
        return result0*2. + result1;
    }
    else {
        VertexInput inputClosedAbove(i0_vertex_right, 0, 0, vp, v, i_in, pick_channel<version>());
#if VECTORIZED_INTEGRATION == 1 and KELDYSH_FORMALISM
        using valuetype_fetch = Eigen::Matrix<Q,4*return_type::ColsAtCompileTime,1>;
        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> result = vertex.template value_symmetry_expanded<pick_spin<version,0>(),pick_channel<version>(),valuetype_fetch>(inputClosedAbove);
        result.resize(4,return_type::ColsAtCompileTime);
#else
        const buffertype_vertex result = vertex.template value_symmetry_expanded<pick_spin<version,0>(),pick_channel<version>(),buffertype_vertex>(inputClosedAbove);
#endif
        return result;
    }
}



#endif //KELDYSH_MFRG_INTEGRANDSE_HPP
