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
* @param[in]  vec_real_offset offsets of the mesh.
*/
void CriterionManager::InitialAllGeometrySerial(const DefInt dims,
    const DefReal reference_dx, std::vector<DefReal> vec_real_offset) {
    numb_of_geometry_ = static_cast<DefInt>(vec_ptr_geometries_.size());
    for (auto i_geo = 0; i_geo < numb_of_geometry_; ++i_geo) {
       // assign offset distance
       bool bool_dims = vec_ptr_geometries_.at(i_geo)->SetOffset(vec_real_offset);
        if (!bool_dims) {
            std::string msg = "dimension of center for geometry " + std::to_string(i_geo) + "is "
                + std::to_string(vec_real_offset.size()) + " but it should be " + std::to_string(dims)
                + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__);
                LogManager::LogError(msg);
        }

        // initial geometry vertex
        vec_ptr_geometries_.at(i_geo)->SetIndex();
        vec_ptr_geometries_.at(i_geo)->InitialGeometry(reference_dx);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
