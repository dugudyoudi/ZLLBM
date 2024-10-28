//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape_2d.cpp
* @author Zhengliang Liu
* @brief functions to initialize and update 2d geometry shapes
* @date  2022-8-5
*/
#include <limits>
#include "./auxiliary_inline_func.h"
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
        LogManager::LogError("pointer to geometry infomation instance is nullptr");
    }
    const DefReal dx = dx_background/TwoPowerN(ptr_geo_info_->GetLevel());
    const std::array<DefReal, 3> geo_center = ptr_geo_info_->GetGeometryCenter(),
        offset = ptr_geo_info_->GetOffset();
    ptr_geo_info_->SetFloodFillOrigin(geo_center);
    DefSizet num_points = DefSizet(2*kPi*radius_ / dx + kEps) + 1;
    ptr_geo_info_->vec_vertices_.resize(num_points);
    DefReal i_real;
    for (DefSizet i = 0; i < num_points; ++i) {
        ptr_geo_info_->vec_vertices_.at(i) = ptr_geo_info_->GeoVertexCreator();
        i_real = (static_cast<DefReal>(i) / static_cast<DefReal>(num_points));
        ptr_geo_info_->vec_vertices_.at(i)->coordinate_[kXIndex] =
            radius_ * cos(2.f * kPi * i_real)
            + geo_center[kXIndex] + offset[kXIndex];
        ptr_geo_info_->vec_vertices_.at(i)->coordinate_[kYIndex] =
            radius_ * sin(2.f * kPi * i_real)
            + geo_center[kYIndex] + offset[kYIndex];
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
    const std::array<DefReal, 3> geo_center = ptr_geo_info_->GetGeometryCenter(),
        offset = ptr_geo_info_->GetOffset();
    ptr_geo_info_->SetFloodFillOrigin(geo_center);
    DefReal length = std::sqrt(Square(end_point_.at(kXIndex) - start_point_.at(kXIndex))
        + Square(end_point_.at(kYIndex) - start_point_.at(kYIndex)));
    DefSizet num_points = DefSizet(length / dx + kEps) + 1;
    ptr_geo_info_->vec_vertices_.resize(num_points);
    DefReal arc = length / num_points,
        sin_theta = (end_point_.at(kYIndex) - start_point_.at(kYIndex)) / length,
        cos_theta = (end_point_.at(kXIndex) - start_point_.at(kXIndex)) / length;
    for (DefSizet i = 0; i < num_points; ++i) {
        ptr_geo_info_->vec_vertices_.at(i) = ptr_geo_info_->GeoVertexCreator();
        ptr_geo_info_->vec_vertices_.at(i)->coordinate_[kXIndex] =
            start_point_.at(kXIndex) + (i + 0.5) * cos_theta * arc + offset[kXIndex];
        ptr_geo_info_->vec_vertices_.at(i)->coordinate_[kYIndex] =
            start_point_.at(kYIndex) + (i + 0.5) * sin_theta * arc + offset[kYIndex];
    }
}
void GeoShapeDefaultLine2D::UpdateShape(const DefReal sum_t) {
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
