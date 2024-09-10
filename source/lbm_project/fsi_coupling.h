//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file fsi_coupling.h
* @author Zhengliang Liu
* @brief define classes to manage fluid structure coupling.
* @date  2024-2-03
*/
#ifndef SOURCE_LBM_PROJECT_FSI_COUPLING_H_
#define SOURCE_LBM_PROJECT_LBM_FSI_COUPLING_H_
#include <vector>
#include <array>
#include <map>
#include "../defs_libs.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace lbmproject {
class GridInfoLbmInteface;
class FsiCoupling {
    int WeakBoundaryCoupling(const DefInt i_level,
        GridInfoLbmInteface* const ptr_grid_info,
        amrproject::GeometryInfoInterface* const ptr_geo_info);
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_FSI_COUPLING_H_
