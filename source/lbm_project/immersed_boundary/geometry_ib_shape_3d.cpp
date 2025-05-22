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
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to update vertices information (area) a 3D quadrilateral
* @param[in] sum_t current time
*/
void GeoShapeIBQuadrilateral3D::UpdateShape(const DefReal sum_t) {
    DefReal edge1_length = std::sqrt(Square(start_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(start_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    DefReal edge2_length = std::sqrt(Square(diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));

    DefReal edge1_arc = edge1_length / edge1_num_points_,
        edge2_arc = edge2_length / edge2_num_points_;
    GeometryInfoImmersedBoundary* ptr_ib_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_ib_geo_info = dynamic_cast<GeometryInfoImmersedBoundary*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    for (auto& iter_vertex : ptr_ib_geo_info->map_vertices_info_) {
        GeometryVertexImmersedBoundary& vertex_ib =
            dynamic_cast<GeometryVertexImmersedBoundary&>(*iter_vertex.second.get());
        vertex_ib.area_ = edge1_arc*edge2_arc;
    }
}
/**
* @brief   function to update vertices information (area) a 3D sphere
* @param[in] sum_t current time
*/
void GeoShapeIBSphere3D::UpdateShape(const DefReal sum_t) {
    // it is assumed the vertices are evenly distributed on the sphere
    GeometryInfoImmersedBoundary* ptr_ib_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_ib_geo_info = dynamic_cast<GeometryInfoImmersedBoundary*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    DefReal area = 4.*kPi*Square(radius_)/ptr_ib_geo_info->vec_vertices_.size();

    for (auto& iter_vertex : ptr_ib_geo_info->map_vertices_info_) {
        GeometryVertexImmersedBoundary& vertex_ib =
            dynamic_cast<GeometryVertexImmersedBoundary&>(*iter_vertex.second.get());
        vertex_ib.area_ = area;
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace lbmproject
}  // end namespace rootproject
