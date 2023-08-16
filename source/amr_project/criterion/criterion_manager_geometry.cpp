//  Copyright (c) 2021 - 2023, Zhengliang Liu
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
* @param[in]  dims dimension of the grid
* @param[in]  vec_real_offset offsets of the mesh
*/
void CriterionManager::InitialAllGeometrySerial(
  const DefAmrIndexUint dims, std::vector<DefReal> vec_real_offset) {
    numb_of_geometry_ = static_cast<DefAmrIndexUint>(vec_ptr_geometries_.size());
    for (auto i_geo = 0; i_geo < numb_of_geometry_; ++i_geo) {
       // assign offset distance
       bool bool_dims = vec_ptr_geometries_.at(i_geo)->SetCenter(vec_real_offset);
       if (!bool_dims) {
        std::string msg = "dimension of center for geometry " + std::to_string(i_geo) + "is "
         + std::to_string(vec_real_offset.size()) + " but it should be " + std::to_string(dims);
         LogManager::LogError(msg);
       }

       // initial geometry vertex
       vec_ptr_geometries_.at(i_geo)->SetIndex();
    //    int initial_msg = vec_ptr_geometries_.at(i_geo)->InitialGeometry(dims, vec_ptr_geometries_.at(i_geo));
    //    if (initial_msg) {
    //        std::string msg = "type of GeometryShape is undefined for geometry " + std::to_string(i_geo);
    //        LogError(msg);
    //    }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
