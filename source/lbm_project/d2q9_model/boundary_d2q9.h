//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9.cpp
* @author Zhengliang Liu
* @brief define class used for LBM D2Q9 model.
* @date  2023-9-30
*/
#ifndef SOURCE_LBM_PROJECT_D2Q9_MODEL_BOUNDARY_D2Q9_H_
#define SOURCE_LBM_PROJECT_D2Q9_MODEL_BOUNDARY_D2Q9_H_
#include <array>
#include <memory>
#include "boundary_conditions/lbm_boundary_conditions.h"
namespace rootproject {
namespace lbmproject {
class GridInfoLbmInteface;
class BoundaryBounceBackD2Q9 : public BoundaryBounceBack2D {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const pr_grid_info) const override;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_D2Q9_MODEL_BOUNDARY_D2Q9_H_
