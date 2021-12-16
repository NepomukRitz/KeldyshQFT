//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP

#include "../parameters/master_parameters.hpp"
#include "../utilities/math_utils.hpp"
#include "Keldysh_symmetries.hpp"

/**
 * specifies the diagrammatic contribution ( + modifications by prefactor, complex conjugation)
 */
struct IndicesSymmetryTransformations{
    int iK;
    double prefactor = 1.; // fermionic sign factor; comes in effect for T1, T2 (and sometimes Tc)
    bool conjugate = false;
    bool asymmetry_transform = false;
    int spin;
    int iw;
    double w, v1, v2; int i_in;
    K_class kClass_aim;     // we only distinguish kClass==k3 from kClass!=k3 --> important for interpolation in K3
    char channel;
    char channel_bubble;
    char channel_parametrization = channel; // W.r.t. which channel is the vertex parametrized? Used for the Hubbard model.

    IndicesSymmetryTransformations(int iK_in, int spin_in, double w_in, double v1_in, double v2_in, int i_in_in, char channel_in, K_class k_in, int iw_in, char channel_bubble_in)
            : iK(iK_in), spin(spin_in), w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in), channel(channel_in), kClass_aim(k_in), iw(iw_in), channel_bubble(channel_bubble_in)
    {}

    IndicesSymmetryTransformations(VertexInput input, char channel_in)
            : iK(input.iK), spin(input.spin), w(input.w), v1(input.v1), v2(input.v2), i_in(input.i_in), channel(channel_in), kClass_aim(input.kClass_aim), iw(input.iw_r), channel_bubble(input.channel)
    {}
};

void switch_channel(IndicesSymmetryTransformations& indices);
void T1 (IndicesSymmetryTransformations& indices);
void T2 (IndicesSymmetryTransformations& indices);
void T3 (IndicesSymmetryTransformations& indices);
void TC (IndicesSymmetryTransformations& indices);
void Tph (IndicesSymmetryTransformations& indices);
void TR (IndicesSymmetryTransformations& indices);
void Ti (IndicesSymmetryTransformations& indices, int i);



#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP