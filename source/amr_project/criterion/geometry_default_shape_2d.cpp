//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape_2d.cpp
* @author Zhengliang Liu
* @brief functions to initialize and update 2d geometry shapes
* @date  2022-8-5
*/
#include <limits>
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to generate a 2d circle
* @param[in] dx_background reference spatial step (background)
*/
void GeoShapeDefaultCircle2D::InitialShape(const DefReal dx_background) {
    if (ptr_geo_info_ == nullptr) {
        LogManager::LogError("pointer to geometry infomation instance is nullptr in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    GeometryInfo2DInterface* ptr_geo = dynamic_cast<GeometryInfo2DInterface*>(ptr_geo_info_);
    const DefReal dx = dx_background/TwoPowerN(ptr_geo->i_level_);
    ptr_geo->flood_fill_origin_ = ptr_geo->geometry_center_;
    DefSizet num_points_ = DefSizet(2*kPi*radius_ / dx + kEps) + 1;
    ptr_geo->coordinate_origin_ = std::vector<GeometryCoordinate2D>(num_points_);
    DefReal i_real;
    for (DefSizet i = 0; i < num_points_; ++i) {
        i_real = (static_cast<DefReal>(i) / static_cast<DefReal>(num_points_));
        ptr_geo->coordinate_origin_.at(i).coordinate[kXIndex] =
            radius_ * cos(2.f * kPi * i_real)
            + ptr_geo->geometry_center_[kXIndex] + ptr_geo->k0RealMin_[kXIndex];
        ptr_geo->coordinate_origin_.at(i).coordinate[kYIndex] =
            radius_ * sin(2.f * kPi * i_real)
            + ptr_geo->geometry_center_[kYIndex] + ptr_geo->k0RealMin_[kYIndex];
    }
}
void GeoShapeDefaultCircle2D::UpdateShape(const DefReal sum_t) {
}
/**
* @brief   function to generate a 2d circle
* @param[in] dx reference spatial step
*/
void GeoShapeDefaultLine2D::InitialShape(const DefReal dx) {
    if (ptr_geo_info_ == nullptr) {
        LogManager::LogError("pointer to geometry infomation instance is nullptr in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    GeometryInfo2DInterface* ptr_geo = dynamic_cast<GeometryInfo2DInterface*>(ptr_geo_info_);
    ptr_geo->flood_fill_origin_ = ptr_geo->geometry_center_;
    DefReal length = std::sqrt(Square(end_point_.at(kXIndex) - start_point_.at(kXIndex))
        + Square(end_point_.at(kYIndex) - start_point_.at(kYIndex)));
    num_points_ = DefSizet(length / dx + kEps) + 1;
    ptr_geo->coordinate_origin_ = std::vector<GeometryCoordinate2D>(num_points_);
    DefReal arc = length / num_points_,
        sin_theta = (end_point_.at(kYIndex) - start_point_.at(kYIndex)) / length,
        cos_theta = (end_point_.at(kXIndex) - start_point_.at(kXIndex)) / length;
    for (DefSizet i = 0; i < num_points_; ++i) {
        ptr_geo->coordinate_origin_.at(i).coordinate[kXIndex] =
            start_point_.at(kXIndex) + (i + 0.5) * cos_theta * arc + ptr_geo->k0RealMin_[kXIndex];
        ptr_geo->coordinate_origin_.at(i).coordinate[kYIndex] =
            start_point_.at(kYIndex) + (i + 0.5) * sin_theta * arc + ptr_geo->k0RealMin_[kYIndex];
    }
}
void GeoShapeDefaultLine2D::UpdateShape(const DefReal sum_t) {
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
