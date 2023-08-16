//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @file default_geo_shape_2d.cpp
* @author Zhengliang Liu
* @brief functions to generate default geomtries
* @date  2022-8-5
*/
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to generate a 2D circle
*/
void DefaultGeoManager::circle_initial(
    GeometryInfo2DInterface* const ptr_geo) const {
    DefReal radius = 0.5;
    DefSizet num_points = 400;
    ptr_geo->flood_fill_origin_ = ptr_geo->geometry_center_;
    ptr_geo->coordinate_origin_ = std::vector<GeometryCoordinate2D>(num_points);
    for (DefSizet i = 0; i < num_points; ++i) {
        ptr_geo->coordinate_origin_.at(i).coordinate[kXIndex] = radius
            * cos(2.f * kPi * (static_cast<DefReal>(i) / static_cast<DefReal>(num_points)))
            + ptr_geo->geometry_center_[kXIndex] + ptr_geo->k0RealOffset_[kXIndex];
        ptr_geo->coordinate_origin_.at(i).coordinate[kYIndex] = radius
            * sin(2.f * kPi * (static_cast<DefReal>(i) / static_cast<DefReal>(num_points)))
            + ptr_geo->geometry_center_[kYIndex] + ptr_geo->k0RealOffset_[kYIndex];
    }
}
void DefaultGeoManager::circle_update(DefReal sum_t,
    GeometryInfo2DInterface* const ptr_geo) const {
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
