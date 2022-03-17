//
// Created by Sa.Aguirre on 7/16/20.
//

#ifndef KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP
#define KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP

#include "../parameters/master_parameters.hpp"
#include "../utilities/math_utils.hpp"
#include "Keldysh_symmetries.hpp"
#include "../data_structures.hpp"


void switch_channel(IndicesSymmetryTransformations& indices);
void T1 (IndicesSymmetryTransformations& indices);
void T2 (IndicesSymmetryTransformations& indices);
void T3 (IndicesSymmetryTransformations& indices);
void TC (IndicesSymmetryTransformations& indices);
void Tph (IndicesSymmetryTransformations& indices);
void TR (IndicesSymmetryTransformations& indices);
void Ti (IndicesSymmetryTransformations& indices, int i);



#endif //KELDYSH_MFRG_SYMMETRY_TRANSFORMATIONS_HPP