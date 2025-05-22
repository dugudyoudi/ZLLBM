//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
#include "io/input_parser.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to generate a 2d circle.
* @param[in] dx_background reference spatial step (background).
*/
void GeoShapeDefaultCircle2D::InitialShape(const DefReal dx_background) {
    GeometryInfoInterface* ptr_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_geo_info = ptr.get();
    } else {
        LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const DefReal dx = dx_background/TwoPowerN(ptr_geo_info->GetLevel());
    const std::array<DefReal, 3> offset = ptr_geo_info->GetOffset();
    ptr_geo_info->SetGeometryCenter({center_[kXIndex]+offset[kXIndex], center_[kYIndex]+offset[kYIndex], 0});
    ptr_geo_info->SetFloodFillOrigin({center_[kXIndex]+offset[kXIndex], center_[kYIndex]+offset[kYIndex], 0});
    DefSizet num_points = DefSizet(2*kPi*radius_ / dx + kEps) + 1;
    ptr_geo_info->vec_vertices_.resize(num_points);
    DefReal i_real;
    for (DefSizet i = 0; i < num_points; ++i) {
        ptr_geo_info->vec_vertices_.at(i) = ptr_geo_info->GeoIndexVertexCreator();
        i_real = (static_cast<DefReal>(i) / static_cast<DefReal>(num_points));
        ptr_geo_info->vec_vertices_.at(i)->coordinate_[kXIndex] =
            radius_ * cos(2. * kPi * i_real) + center_[kXIndex] + offset[kXIndex];
        ptr_geo_info->vec_vertices_.at(i)->coordinate_[kYIndex] =
            radius_ * sin(2. * kPi * i_real) + center_[kYIndex] + offset[kYIndex];
    }
}
/**
* @brief   function to read and set parameters for a 2d circle.
* @param[in] shape_parameters map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultCircle2D::ReadAndSetGeoShapeParameters(
    std::map<std::string, ParserData>* const ptr_shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("cylinder.center", ptr_shape_parameters, &center_)) {
        LogManager::LogWarning("center is not found for 2d circle");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("cylinder.radius", ptr_shape_parameters, &radius_)) {
        LogManager::LogWarning("radius is not found for 2d circle");
        parameter_exist = false;
    }
    return parameter_exist;
}
/**
* @brief   function to generate a 2d line
* @param[in] dx reference spatial step
*/
void GeoShapeDefaultLine2D::InitialShape(const DefReal dx) {
    GeometryInfoInterface* ptr_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_geo_info = ptr.get();
    } else {
        LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const std::array<DefReal, 3> offset = ptr_geo_info->GetOffset();
    std::array<DefReal, 2> center = {0.5*(start_point_.at(kXIndex) + end_point_.at(kXIndex)),
        0.5*(start_point_.at(kYIndex) + end_point_.at(kYIndex))};
    ptr_geo_info->SetGeometryCenter({center[kXIndex]+offset[kXIndex], center[kYIndex]+offset[kYIndex], 0});
    ptr_geo_info->SetFloodFillOrigin({center[kXIndex]+offset[kXIndex], center[kYIndex]+offset[kYIndex], 0});
    DefReal length = std::sqrt(Square(end_point_.at(kXIndex) - start_point_.at(kXIndex))
        + Square(end_point_.at(kYIndex) - start_point_.at(kYIndex)));
    DefSizet num_points = DefSizet(length / dx + kEps);
    ptr_geo_info->vec_vertices_.resize(num_points);
    DefReal arc = length / num_points,
        sin_theta = (end_point_.at(kYIndex) - start_point_.at(kYIndex)) / length,
        cos_theta = (end_point_.at(kXIndex) - start_point_.at(kXIndex)) / length;
    for (DefSizet i = 0; i < num_points; ++i) {
        ptr_geo_info->vec_vertices_.at(i) = ptr_geo_info->GeoIndexVertexCreator();
        ptr_geo_info->vec_vertices_.at(i)->coordinate_[kXIndex] =
            start_point_.at(kXIndex) + (i + 0.5) * cos_theta * arc + offset[kXIndex];
        ptr_geo_info->vec_vertices_.at(i)->coordinate_[kYIndex] =
            start_point_.at(kYIndex) + (i + 0.5) * sin_theta * arc + offset[kYIndex];
    }
}
/**
* @brief   function to read and set parameters for a 2d line.
* @param[in] ptr_shape_parameters pointer to map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultLine2D::ReadAndSetGeoShapeParameters(
    std::map<std::string, ParserData>* const ptr_shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("line.start_point", ptr_shape_parameters, &start_point_)) {
        LogManager::LogWarning("start point is not found for 2d line");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("line.end_point", ptr_shape_parameters, &end_point_)) {
        LogManager::LogWarning("end point is not found for 2d line");
        parameter_exist = false;
    }
    return parameter_exist;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
