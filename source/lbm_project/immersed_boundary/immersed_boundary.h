//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.h
* @author Zhengliang Liu
* @brief define classes to manage immersed boundary method.
* @date  2024-2-03
*/
#ifndef ROOTPROJECT_SOURCE_LBM_IMMERSED_BOUNDARY_H_
#define ROOTPROJECT_SOURCE_LBM_IMMERSED_BOUNDARY_H_
#include <vector>
#include <array>
#include <map>
#include "fsi_coupling.h"
namespace rootproject {
namespace lbmproject {
class FsiImmersedBoundary : public FsiCoupling {
    DefReal StencilDisOne(DefReal dis);
    DefReal StencilDisTwo(DefReal dis);
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_IMMERSED_BOUNDARY_H_
