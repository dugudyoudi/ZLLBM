//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage changes of geometries
* @date  2022-5-23
* @note .
*/
#include "grid/grid_manager.h"
#include "io/log_write.h"
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to initialize geometries
* @param[in]  dims dimension of the grid.
* @param[in]  reference_dx reference spatial step.
* @param[in]  real_offset offsets of the mesh.
*/
void CriterionManager::InitialAllGeometrySerial(const DefInt dims,
    const DefReal reference_dx, const std::array<DefReal, 3>& real_offset) {
    numb_of_geometry_ = static_cast<DefInt>(vec_ptr_geometries_.size());
    for (auto i_geo = 0; i_geo < numb_of_geometry_; ++i_geo) {
        // assign offset distance
        vec_ptr_geometries_.at(i_geo)->SetOffset(real_offset);

        // initial geometry vertex
        vec_ptr_geometries_.at(i_geo)->InitialGeometry(reference_dx);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
