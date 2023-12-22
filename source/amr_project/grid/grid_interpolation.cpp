//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_interpolation.cpp
* @author Zhengliang Liu
* @brief functions used for node interpolation.
* @date  2023-12-23
*/
#include <vector>
#include <map>
#include "grid/grid_info_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
const GridInfoInterface::LagrangianCoeff& GridInfoInterface::CalculateLagrangianInterpCoeff(
    const DefAmrIndexLUint interp_half_length) {
    if (lagrangian_coefficients_.find(interp_half_length) == lagrangian_coefficients_.end()) {
        DefAmrIndexLUint offset = interp_half_length - 1;
        lagrangian_coefficients_.insert({interp_half_length, LagrangianCoeff()});
        DefAmrIndexLUint num_coff = 2 * interp_half_length;
        lagrangian_coefficients_.at(interp_half_length).coeff0.assign(num_coff, 0.);
        lagrangian_coefficients_.at(interp_half_length).coeff0.shrink_to_fit();
        lagrangian_coefficients_.at(interp_half_length).coeff1.assign(num_coff, 0.);
        lagrangian_coefficients_.at(interp_half_length).coeff1.shrink_to_fit();
        DefReal x_tmp, y_tmp;
        for (DefAmrIndexLUint iy = 0; iy < num_coff; ++iy) {
            lagrangian_coefficients_.at(interp_half_length).coeff1.at(iy) = 1.;
            y_tmp = static_cast<DefReal>(iy);
            for (DefAmrIndexLUint ix = 0; ix < num_coff; ++ix) {
                if (iy != ix) {
                    x_tmp = static_cast<DefReal>(ix);
                    lagrangian_coefficients_.at(interp_half_length).coeff1.at(iy) *=
                        (0.5 + offset - x_tmp)/(y_tmp - x_tmp);
                }
            }
        }
        lagrangian_coefficients_.at(interp_half_length).coeff0.at(offset) = 1.;
    }
    return lagrangian_coefficients_.at(interp_half_length);
}

}  // end namespace amrproject
}  // end namespace rootproject
