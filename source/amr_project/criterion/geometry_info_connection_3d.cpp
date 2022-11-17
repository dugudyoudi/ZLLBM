//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_stl.cpp
* @author Zhengliang Liu
* @brief functions for geometries reprensented by STL vertexs
* @date  2022-5-25
* @note .
*/
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_connection.h"
#include "criterion/criterion_manager.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
namespace criterion {
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
* @note
*/
void GeometryInfoConnection3D::SetIndex() {
    if (bool_vec_velocities_) {
        k0UxIndex_ = 0; k0UyIndex_ = k0UxIndex_ + 1;
        k0UzIndex_ = k0UyIndex_ + 1;
        k0NumRealForEachvertex_ = k0UzIndex_;
    }
    if (bool_vec_forces_) {
        k0FxIndex_ = k0UzIndex_ + 1; k0FyIndex_ = k0FxIndex_ + 1;
        k0FzIndex_ = k0FyIndex_ + 1;
        k0NumRealForEachvertex_ = k0FzIndex_;
    }
    vertex_instance_.coordinates = { 0., 0., 0.};
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachvertex_);
}
/**
* @brief   function to initialize status of geometries
* @return  0 successful, 1 geometry is undefined
*/
int GeometryInfoConnection3D::InitialGeometry(
    std::shared_ptr<GeometryInfo3DInterface> ptr_geo) {
    DefUint dims = 3;
    SetIndex();

    return 0;
}
int GeometryInfoConnection3D::UpdateGeometry(
    std::shared_ptr<GeometryInfo3DInterface> ptr_geo) {
    return 0;
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
*/
void GeometryInfoConnection3D::InitialCoordinateGivenLevel() {
    vertex_given_level_.push_back({});
    connection_vertex_given_level_.push_back({});
    DefSizet i_vertex = 0;
    GeometryConnectionCoordinate vertex_temp;
    vertex_temp.coordinates = { 0., 0., 0. };
    for (const auto& iter_vertex : coordinate_origin_) {
        vertex_temp.coordinates.at(0) = iter_vertex.coordinate.at(0);
        vertex_temp.coordinates.at(1) = iter_vertex.coordinate.at(1);
        vertex_temp.coordinates.at(2) = iter_vertex.coordinate.at(2);
        vertex_given_level_.at(0).vec_vertex_cooridinate
            .push_back(vertex_temp);
        connection_vertex_given_level_.at(0).insert({ 0, i_vertex });
        ++i_vertex;
    }
}
/**
* @brief function to compute distance between two vertices (3D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @return  distance.
*/
DefReal GeometryInfoConnection3D::ComputeDistanceFromCoordinates(
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
    DefReal z_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kZIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_cooridinate.at(vertex1.second).coordinates.at(kZIndex));
    return sqrt(x_dis * x_dis + y_dis * y_dis + z_dis * z_dis);
}
/**
* @brief function to compute coordinates of the mid point of two vertices (3D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @param[out]  ptr_coordinates mid point coordinates.
*/
void GeometryInfoConnection3D::ComputeMidCoordinates(
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
    ptr_coordinates->at(kZIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_cooridinate.at(vertex0.second).coordinates.at(kZIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_cooridinate
        .at(vertex1.second).coordinates.at(kZIndex)) / 2.;
}
void GeometryInfoConnection3D::DecomposeNHigerLevel(const DefSizet i_level_grid,
    const DefReal decompose_length,
    const std::unordered_map<DefSizet, bool>& map_indices_base,
    std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain){

};
void GeometryInfoConnection3D::FindTrackingNodeNearGeo(
    const grid::SFBitsetAux3D& sfbitset_aux_3d,
    std::shared_ptr<grid::GridInfoInterface> ptr_grid_info) {

}
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_3D_FUNCTIONS