//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_stl.cpp
* @author Zhengliang Liu
* @brief functions for geometries reprensented by STL vertexs
* @date  2022-5-25
* @note .
*/
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
#include <algorithm>
//#include <stack>
#include <queue>
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/criterion_manager.h"
#include "criterion/geometry_info_connection.h"
#include "grid/grid_info_interface.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
*/
void GeometryInfoConnection2D::SetIndex() {
    if (bool_vec_velocities_) {
        k0UxIndex_ = 0; k0UyIndex_ = k0UxIndex_ + 1;
        k0NumRealForEachvertex_ = k0UyIndex_;
    }
    if (bool_vec_forces_) {
        k0FxIndex_ = k0UyIndex_ + 1; k0FyIndex_ = k0FxIndex_ + 1;
        k0NumRealForEachvertex_ = k0FyIndex_;
    }
    vertex_instance_.coordinates = { 0., 0. };
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachvertex_);
}
/**
* @brief   function to initialize status of geometries
* @return  0 successful, 1 geometry is undefined
* @note
*/
int GeometryInfoConnection2D::InitialGeometry(
    std::shared_ptr<GeometryInfo2DInterface> ptr_geo) {
    DefUint dims = 2;
    SetIndex();
    SetupConnectionParameters(geometry_cell_type_);
    switch (k0DefaultGeoShapeType_) {
    case DefaultGeoShapeType::kCircle:
        DefaultGeoShapeManager::GetInstance()
            ->geo_circle_initial(dims, ptr_geo);
        return 0;
    default:
        return 1;
    }
}
int GeometryInfoConnection2D::UpdateGeometry(
    std::shared_ptr<GeometryInfo2DInterface> ptr_geo) {
    return 0;
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
*/
void GeometryInfoConnection2D::InitialCoordinateGivenLevel() {
    vertex_given_level_.push_back({});
    connection_vertex_given_level_.push_back({});
    DefSizet i_vertex = 0;
    GeometryConnectionCoordinate vertex_temp;
    vertex_temp.coordinates = { 0., 0.};
    for (const auto& iter_vertex : coordinate_origin_) {
        vertex_temp.coordinates.at(0) = iter_vertex.coordinate.at(0);
        vertex_temp.coordinates.at(1) = iter_vertex.coordinate.at(1);
        vertex_given_level_.at(0).vec_vertex_cooridinate
            .push_back(vertex_temp);
        connection_vertex_given_level_.at(0).insert({ 0, i_vertex });
        vertex_given_level_.at(0).vec_vertex_cooridinate.at(i_vertex)
            .map_linked_vertices_level.insert({});
        ++i_vertex;
    }
}
/**
* @brief function to compute distance between two vertices (2D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @return  distance.
*/
DefReal GeometryInfoConnection2D::ComputeDistanceFromCoordinates(
    const std::pair<DefSizet, DefSizet>& vertex0,
    const std::pair<DefSizet, DefSizet>& vertex1) {
    DefReal x_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kXIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_cooridinate.at(vertex1.second).coordinates.at(kXIndex));
    DefReal y_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kYIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_cooridinate.at(vertex1.second).coordinates.at(kYIndex));
    return sqrt(x_dis * x_dis + y_dis * y_dis);
}
/**
* @brief function to compute coordinates of the mid point of two vertices (2D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @param[out]  ptr_coordinates mid point coordinates.
*/
void GeometryInfoConnection2D::ComputeMidCoordinates(
    const std::pair<DefSizet, DefSizet>& vertex0,
    const std::pair<DefSizet, DefSizet>& vertex1,
    std::vector<DefReal>* const ptr_coordinates) {
    ptr_coordinates->at(kXIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kXIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_cooridinate
        .at(vertex1.second).coordinates.at(kXIndex)) / 2.;
    ptr_coordinates->at(kYIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kYIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_cooridinate
        .at(vertex1.second).coordinates.at(kYIndex)) / 2.;
}

void GeometryInfoConnection2D::DecomposeNHigerLevel(const DefSizet i_level_grid,
    const DefReal decompose_length,
    const std::unordered_map<DefSizet, bool>& map_indices_base,
    std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) {

}
void GeometryInfoConnection2D::FindTrackingNodeNearGeo(
    const SFBitsetAux2D& sfbitset_aux_2d,
    std::shared_ptr<GridInfoInterface> ptr_grid_info) const {
    if (ptr_grid_info->map_ptr_tracking_grid_info_
        .find({ ECriteriolType::kGeometry, i_geo_ })
        == ptr_grid_info->map_ptr_tracking_grid_info_.end()) {
        ptr_grid_info->map_ptr_tracking_grid_info_.insert(
            {{ ECriteriolType::kGeometry, i_geo_ },
            ptr_tracking_grid_info_creator_->CreateTrackingGridInfo()});
    }
    TrackingGridInfoInterface* ptr_tracking_grid = ptr_grid_info->
        map_ptr_tracking_grid_info_.at(
            { ECriteriolType::kGeometry, i_geo_ }).get();
    DefSizet level_diff = ptr_grid_info->i_level_ - i_level_;
    std::array<DefReal, 2> coordinates;
    for (const auto& iter : connection_vertex_given_level_.at(level_diff)) {
        coordinates[kXIndex] = vertex_given_level_.at(iter.first)
            .vec_vertex_cooridinate.at(iter.second).coordinates[kXIndex];
        //sfbitset_aux_2d.SFBitsetEncoding(coordinates);
           
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
