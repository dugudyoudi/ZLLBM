//  Copyright (c) 2022, Zhengliang Liu
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
* @brief   function to create instance based on the type of geometry vertex
* @param[in]  geo_data_type date type of the geometry
*/
void CriterionManager::CreateGeometryInfo(const DefUint i_level,
    const DefUint geo_data_type) {
    ++numb_of_geometry_;
    vec_ptr_geometries_.back()->node_type_ = geo_data_type;
    // set default i_level_ in GeometryInfo as GridManager::k0MaxLevel_
    vec_ptr_geometries_.back()->i_level_ = i_level;
}
/**
* @brief   function to initialize geometries
*/
void CriterionManager::InitialGeometrySerial(
    std::vector<DefReal> vec_real_offset) {
    numb_of_geometry_ = vec_ptr_geometries_.size();
    //for (auto i_geo = 0; i_geo < numb_of_geometry_; ++i_geo) {
    //    // check if dimension of the given geometry center is correct
    //    if (k0GeoDims_ != vec_ptr_geometries_.at(i_geo)->geometry_center_.size()) {
    //        std::string msg = "Dimension of center coordinate for geometry "
    //            + std::to_string(i_geo)
    //            + " is different from dimension (k0GeoDims_), reset center at ("
    //            + std::to_string
    //            (vec_ptr_geometries_.at(i_geo)->geometry_center_[0]) + ", "
    //            + std::to_string
    //            (vec_ptr_geometries_.at(i_geo)->geometry_center_[1]);
    //        if (k0GeoDims_ == 2) {
    //            msg += ").";
    //        } else {
    //            msg += std::to_string
    //            (vec_ptr_geometries_.at(i_geo)->geometry_center_[2]) + ").";
    //        }
    //        io::LogWarning(msg);
    //        vec_ptr_geometries_.at(i_geo)->geometry_center_.resize(k0GeoDims_);
    //        vec_ptr_geometries_.at(i_geo)->geometry_center_.shrink_to_fit();
    //    }

    //    // assign offset distance
    //    vec_ptr_geometries_.at(i_geo)->k0RealOffset_[kXIndex]
    //        = vec_real_offset[kXIndex];
    //    vec_ptr_geometries_.at(i_geo)->k0RealOffset_[kYIndex]
    //        = vec_real_offset[kYIndex];
    //    if (k0GeoDims_ == 3) {
    //        vec_ptr_geometries_.at(i_geo)->k0RealOffset_[kZIndex]
    //            = vec_real_offset[kZIndex];
    //    }

    //    // initial geometry vertex 
    //    vec_ptr_geometries_.at(i_geo)->SetIndex();
    //    int initial_msg = vec_ptr_geometries_.at(i_geo)->InitialGeometry(
    //        k0GeoDims_, vec_ptr_geometries_.at(i_geo));
    //    if (initial_msg) {
    //        std::string msg = "type of GeometryShape is undefined"
    //            " for geometry " + std::to_string(i_geo);
    //        io::LogError(msg);
    //    }
    //}
}
}  // end namespace amrproject
}  // end namespace rootproject
