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

/// possible tests: ---> See bubble integrand

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 */
template <typename Q, vertexType vertType>
class IntegrandSE {
    const char type;
    std::vector<int> components = std::vector<int>(6);

    const GeneralVertex<Q,vertType>& vertex;
    const Propagator<Q>& propagator;

    const int iK;
    const double v;
    const int i_in;
    const int i_spin;

    void set_Keldysh_components_to_be_calculated();

    Q Keldysh_value(double vp) const;
    Q Matsubara_value(double vp) const;

    void evaluate_propagator(Q& Gi, const int iK, const double vp) const;
    void evaluate_propagator(Q& GM, const double vp) const; // Matsubara version

    void evaluate_vertex(Q& factorClosedAbove, Q& factorAClosedBelow,
                         const int iK, const double vp) const; // for symmetrized Keldysh flow
    void evaluate_vertex(Q& factorClosedAbove,
                         const int iK, const double vp) const; // for unsymmetrized Keldysh/Matsubara flow

public:
    IntegrandSE(const char type_in, const GeneralVertex<Q,vertType>& vertex_in, const Propagator<Q>& prop_in,
                const int iK_in, const int i_spin_in, const double v_in, const int i_in_in)
                :type(type_in), vertex(vertex_in), propagator(prop_in), iK(iK_in), i_spin(i_spin_in), v(v_in), i_in(i_in_in){
        if (KELDYSH){set_Keldysh_components_to_be_calculated();}
    }

    auto operator()(double vp) const -> Q;

    void save_integrand() const;
    void save_integrand(const rvec& freqs) const;
    void get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im)  const;

};

template<typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::set_Keldysh_components_to_be_calculated() {
    if(type=='r'){  //Check which kind of contribution is calculated
        components[0]=3;    //Vertex component associated to Retarded propagator
        components[1]=6;    //Vertex component associated to Advanced propagator
        components[2]=7;    //Vertex component associated to Keldysh propagator
        components[3]=3;    //Vertex component associated to Retarded propagator in symmetrized flow
        components[4]=9;    //Vertex component associated to Advanced propagator in symmetrized flow
        components[5]=11;   //Vertex component associated to Keldysh propagator in symmetrized flow
    }
    else {
        components[0]=1;    //Vertex component associated to Retarded propagator
        components[1]=4;    //Vertex component associated to Advanced propagator
        components[2]=5;    //Vertex component associated to Keldysh propagator
        components[3]=2;    //Vertex component associated to Retarded propagator in symmetrized flow
        components[4]=8;    //Vertex component associated to Advanced propagator in symmetrized flow
        components[5]=10;   //Vertex component associated to Keldysh propagator in symmetrized flow
    }
}

template<typename Q, vertexType vertType>
auto IntegrandSE<Q,vertType>::operator()(const double vp) const -> Q {
    if (KELDYSH){return Keldysh_value(vp);}
    else{return Matsubara_value(vp);}
}

template<typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im) const {
    int npoints = freqs.size();
    for (int i=0; i<npoints; ++i) {

        double vpp = freqs[i];


        Q integrand_value;

        integrand_value = (*this)(vpp);

        if (PARTICLE_HOLE_SYMMETRY && (!KELDYSH)){
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
        }
        else{
            integrand_re[i] = myreal(integrand_value);
            integrand_im[i] = myimag(integrand_value);
        }
    }


}

template<typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = 1e5;

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;

        wl = propagator.selfenergy.frequencies.w_lower * 2.;
        wu = propagator.selfenergy.frequencies.w_upper * 2.;

        double vpp = wl + i * (wu - wl) / (npoints - 1);
        freqs[i] = vpp;
    }

    save_integrand(freqs);

}

template<typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::save_integrand(const rvec& freqs) const {
    int npoints = freqs.size();

    rvec integrand_re (npoints);
    rvec integrand_im (npoints);

    get_integrand_vals(freqs, integrand_re, integrand_im);

    std::string filename = data_dir + "integrand_SE";
    filename += //"_i0=" + std::to_string(i0)       /// TODO: add this when Elias interchanged order of integration and Keldysh sum
                //+ "_i2=" + std::to_string(i2)
                + "_v=" + std::to_string(v);
    filename += + ".h5";
    write_h5_rvecs(filename,
                   {"v", "integrand_re", "integrand_im"},
                   {freqs, integrand_re, integrand_im});
}

template<typename Q, vertexType vertType>
Q IntegrandSE<Q,vertType>::Keldysh_value(const double vp) const {
    Q Gi;
    evaluate_propagator(Gi, iK, vp);

    Q factorClosedAbove;
#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
    Q factorClosedBelow;
    evaluate_vertex(factorClosedAbove, factorClosedBelow, iK, vp);
    return (1./2.) * Gi * (factorClosedAbove + factorClosedBelow);
#else
    evaluate_vertex(factorClosedAbove, iK, vp);
    return Gi * factorClosedAbove;
#endif


}

template<typename Q, vertexType vertType>
Q IntegrandSE<Q,vertType>::Matsubara_value(const double vp) const {
    Q GM;
    evaluate_propagator(GM, vp);

    Q factorClosedAbove;
    evaluate_vertex(factorClosedAbove, 0, vp);

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

template <typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::evaluate_propagator(Q &Gi, const int iK, const double vp) const {
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
}

template<typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::evaluate_propagator(Q &GM, const double vp) const {
    GM = propagator.valsmooth(0, vp, i_in);           // Matsubara propagator (full or single scale)
}

template <typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::evaluate_vertex(Q &factorClosedAbove, Q &factorClosedBelow,
                                     const int iK, const double vp) const {
    VertexInput inputClosedAbove (components[iK]  , i_spin, 0., vp, v, i_in, 't');
    VertexInput inputClosedBelow (components[iK+3], i_spin, 0., v, vp, i_in, 't');
    factorClosedAbove = vertex.template value<'t'>(inputClosedAbove);
    factorClosedBelow = vertex.template value<'t'>(inputClosedBelow);
}

template <typename Q, vertexType vertType>
void IntegrandSE<Q,vertType>::evaluate_vertex(Q &factorClosedAbove, const int iK, const double vp) const {
    // "components" are all zero in Matsubara case -> this function also works for Matsubara
    VertexInput inputClosedAbove (components[iK], i_spin, 0, vp, v, i_in, 't');
    factorClosedAbove = vertex.template value<'t'>(inputClosedAbove);
}



#endif //KELDYSH_MFRG_INTEGRANDSE_HPP
