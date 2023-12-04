//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9.cpp
* @author Zhengliang Liu
* @brief define class used for LBM D2Q9 model.
* @date  2023-9-30
*/
#ifndef ROOTPROJECT_SOURCE_LBM_BOUNDARY_D2Q9_H_
#define ROOTPROJECT_SOURCE_LBM_BOUNDARY_D2Q9_H_
#include <array>
#include <memory>
#include "boundary_conditions/lbm_boundary_conditions.h"
namespace rootproject {
namespace lbmproject {
class GridInfoLbmInteface;
class BoundaryBounceBackD2Q9 : public BoundaryBounceBack2D {
 public:
    void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
        GridInfoLbmInteface* const pr_grid_info) const override;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_BOUNDARY_D2Q9_H_
