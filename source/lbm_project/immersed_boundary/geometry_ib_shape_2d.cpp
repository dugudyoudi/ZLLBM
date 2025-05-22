//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_ib_shape_2d.cpp
* @author Zhengliang Liu
* @brief functions to update 2d geometry shapes for immersed boundary method
* @date  2022-8-5
*/
#include "io/log_write.h"
#include "immersed_boundary/geometry_ib_shape.h"
#include "immersed_boundary/immersed_boundary.h"
namespace rootproject {
namespace lbmproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to update vertices information (area) a 2D line
* @param[in] sum_t current time
*/
void GeoShapeIBLine2D::UpdateShape(const DefReal sum_t) {
    const DefReal length = std::sqrt(Square(end_point_.at(kXIndex) - start_point_.at(kXIndex))
        + Square(end_point_.at(kYIndex) - start_point_.at(kYIndex)));
    GeometryInfoImmersedBoundary* ptr_ib_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_ib_geo_info = dynamic_cast<GeometryInfoImmersedBoundary*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const DefSizet num_points = ptr_ib_geo_info->vec_vertices_.size();
    DefReal arc = length / num_points;
    for (auto& iter_vertex : ptr_ib_geo_info->map_vertices_info_) {
        GeometryVertexImmersedBoundary& vertex_ib =
            dynamic_cast<GeometryVertexImmersedBoundary&>(*iter_vertex.second.get());
        vertex_ib.area_ = arc;
    }
}
/**
* @brief   function to update vertices information (area) a 2D circle
* @param[in] sum_t current time
*/
void GeoShapeIBCircle2D::UpdateShape(const DefReal sum_t) {
    GeometryInfoImmersedBoundary* ptr_ib_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_ib_geo_info = dynamic_cast<GeometryInfoImmersedBoundary*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const DefSizet num_points = ptr_ib_geo_info->vec_vertices_.size();
    DefReal arc = 2*kPi*radius_ / num_points;
    for (auto& iter_vertex : ptr_ib_geo_info->map_vertices_info_) {
        GeometryVertexImmersedBoundary& vertex_ib =
            dynamic_cast<GeometryVertexImmersedBoundary&>(*iter_vertex.second.get());
        vertex_ib.area_ = arc;
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}  // end namespace lbmproject
}  // end namespace rootproject
