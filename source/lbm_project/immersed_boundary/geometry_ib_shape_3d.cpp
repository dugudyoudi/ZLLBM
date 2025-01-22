//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
* @brief   function to generate a 2d circle
* @param[in] dx reference spatial step
*/
void GeoShapeIBQuadrilateral3D::UpdateShape(const DefReal sum_t) {
    if (ptr_geo_info_ == nullptr) {
        amrproject::LogManager::LogError("pointer to geometry infomation instance is nullptr in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefReal edge1_length = std::sqrt(Square(start_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(start_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    DefReal edge2_length = std::sqrt(Square(diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));

    DefReal edge1_arc = edge1_length / edge1_num_points_,
        edge2_arc = edge2_length / edge2_num_points_;
    GeometryInfoImmersedBoundary* ptr_ib_geo_info = dynamic_cast<GeometryInfoImmersedBoundary*>(ptr_geo_info_);
    for (auto& iter_vertex : ptr_ib_geo_info->map_vertices_info_) {
        GeometryVertexImmersedBoundary& vertex_ib =
            dynamic_cast<GeometryVertexImmersedBoundary&>(*iter_vertex.second.get());
        vertex_ib.area_ = edge1_arc*edge2_arc;
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}  // end namespace lbmproject
}  // end namespace rootproject
