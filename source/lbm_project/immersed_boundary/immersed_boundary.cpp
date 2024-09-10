//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.cpp
* @author Zhengliang Liu
* @brief functions used to implement immersed boundary method.
* @date  2023-11-6
*/
#include "immersed_boundary/immersed_boundary.h"
namespace rootproject {
namespace lbmproject {
void FsiImmersedBoundary::DirectForcingScheme() {

}

DefReal FsiImmersedBoundary::StencilDisOne(DefReal dist) {
    if (std::fabs(dist) < 1.) {
        return 1. - std::fabs(dist);
    } else {
        return 0.;
    }
}
DefReal FsiImmersedBoundary::StencilDisTwo(DefReal dist) {
    if (dist < 2.) {
        DefReal dist_abs = std::fabs(dist);
        DefReal dist_sq = dist*dist;
        if (dist_abs < 1.) {
            return 0.125 * (3. - 2. * dist_abs + sqrt(1. + 4. * dist_abs - 4. * dist_sq));
        } else {
            return 0.125 * (5. - 2. * dist_abs - sqrt(-7. + 12. * dist_abs - 4. * dist_sq));
        }
    } else {
        return 0.;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
