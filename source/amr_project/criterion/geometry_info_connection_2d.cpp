//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection_2d.cpp
* @author Zhengliang Liu
* @brief functions for 2D geometries with connection information
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
        k0NumRealForEachVertex_ = k0UyIndex_;
    }
    if (bool_vec_forces_) {
        k0FxIndex_ = k0UyIndex_ + 1; k0FyIndex_ = k0FxIndex_ + 1;
        k0NumRealForEachVertex_ = k0FyIndex_;
    }
    vertex_instance_.coordinates = { 0., 0. };
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachVertex_);

    if (bool_vertex_info_stored_for_connection_) {
        vertex_instance_.vertex_info.vec_int.resize(k0NumIntForEachVertex_);
        vertex_instance_.vertex_info.vec_real.resize(k0NumRealForEachVertex_);
    }
}
/**
* @brief   function to initialize status of geometries
* @param[in] dx reference spatial step
* @return  0 successful
*/
int GeometryInfoConnection2D::InitialGeometry(const DefReal dx) {
    int return_status = 0;
    this->SetIndex();
    this->SetupConnectionParameters(this->geometry_cell_type_);
    return_status = this->GeometryInfo2DInterface::InitialGeometry(dx);
    return return_status;
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
* @param[out]  ptr_coordi_min   minimum coordinates of the geometry.
* @param[out]  ptr_coordi_max   maximum coordinates of the geometry.
*/
void GeometryInfoConnection2D::InitialCoordinateGivenLevel(
    std::vector<DefReal>* const ptr_coordi_min,
    std::vector<DefReal>* const ptr_coordi_max) {
    vertex_given_level_.push_back({});
    connection_vertex_given_level_.push_back({});
    DefSizet i_vertex = 0;
    GeometryConnectionCoordinate vertex_tmp;
    vertex_tmp.coordinates = { 0., 0.};
    ptr_coordi_min->resize(2);
    ptr_coordi_max->resize(2);
    ptr_coordi_min->at(kXIndex) =
        coordinate_origin_.at(0).coordinate.at(kXIndex);
    ptr_coordi_min->at(kYIndex) =
        coordinate_origin_.at(0).coordinate.at(kYIndex);
    ptr_coordi_max->at(kXIndex) =
        coordinate_origin_.at(0).coordinate.at(kXIndex);
    ptr_coordi_max->at(kYIndex) =
        coordinate_origin_.at(0).coordinate.at(kYIndex);
    for (const auto& iter_vertex : coordinate_origin_) {
        if (ptr_coordi_min->at(kXIndex) >
            iter_vertex.coordinate.at(kXIndex)) {
            ptr_coordi_min->at(kXIndex) =
                iter_vertex.coordinate.at(kXIndex);
        } else if (ptr_coordi_max->at(kXIndex) <
            iter_vertex.coordinate.at(kXIndex)) {
            ptr_coordi_max->at(kXIndex) =
                iter_vertex.coordinate.at(kXIndex);
        }
        if (ptr_coordi_min->at(kYIndex) >
            iter_vertex.coordinate.at(kYIndex)) {
            ptr_coordi_min->at(kYIndex) =
                iter_vertex.coordinate.at(kYIndex);
        } else if (ptr_coordi_max->at(kYIndex) <
            iter_vertex.coordinate.at(kYIndex)) {
            ptr_coordi_max->at(kYIndex) =
                iter_vertex.coordinate.at(kYIndex);
        }
        vertex_tmp.coordinates.at(kXIndex) = iter_vertex.coordinate.at(kXIndex);
        vertex_tmp.coordinates.at(kYIndex) = iter_vertex.coordinate.at(kYIndex);
        vertex_given_level_.at(0).vec_vertex_coordinate.push_back(vertex_tmp);
        connection_vertex_given_level_.at(0).insert({ 0, i_vertex });
        vertex_given_level_.at(0).vec_vertex_coordinate.at(i_vertex)
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
    const std::pair<DefInt, DefSizet>& vertex0,
    const std::pair<DefInt, DefSizet>& vertex1) {
    DefReal x_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kXIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate.at(vertex1.second).coordinates.at(kXIndex));
    DefReal y_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kYIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate.at(vertex1.second).coordinates.at(kYIndex));
    return sqrt(x_dis * x_dis + y_dis * y_dis);
}
/**
* @brief function to compute coordinates of the mid point of two vertices (2D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @param[out]  ptr_coordinates mid point coordinates.
*/
void GeometryInfoConnection2D::ComputeMidCoordinates(
    const std::pair<DefInt, DefSizet>& vertex0,
    const std::pair<DefInt, DefSizet>& vertex1,
    std::vector<DefReal>* const ptr_coordinates) {
    ptr_coordinates->at(kXIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kXIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate
        .at(vertex1.second).coordinates.at(kXIndex)) / 2.;
    ptr_coordinates->at(kYIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kYIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate
        .at(vertex1.second).coordinates.at(kYIndex)) / 2.;
}

}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
