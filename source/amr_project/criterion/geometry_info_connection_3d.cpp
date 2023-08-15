//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_stl.cpp
* @author Zhengliang Liu
* @brief functions for 3D geometries with connection information
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
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
*/
void GeometryInfoConnection3D::SetIndex() {
    if (bool_vec_velocities_) {
        k0UxIndex_ = 0; k0UyIndex_ = k0UxIndex_ + 1;
        k0UzIndex_ = k0UyIndex_ + 1;
        k0NumRealForEachVertex_ = k0UzIndex_;
    }
    if (bool_vec_forces_) {
        k0FxIndex_ = k0UzIndex_ + 1; k0FyIndex_ = k0FxIndex_ + 1;
        k0FzIndex_ = k0FyIndex_ + 1;
        k0NumRealForEachVertex_ = k0FzIndex_;
    }
    vertex_instance_.coordinates = { 0., 0., 0.};
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachVertex_);

    if (bool_vertex_info_stored_for_connection_) {
        vertex_instance_.vertex_info.vec_int.resize(k0NumIntForEachVertex_);
        vertex_instance_.vertex_info.vec_real.resize(k0NumRealForEachVertex_);
    }
}
/**
* @brief   function to initialize status of geometries
* @return  0 successful, 1 geometry is undefined
*/
int GeometryInfoConnection3D::InitialGeometry(
    const DefReal dx,
    const DefaultGeoShapeType shape_type,
    const DefaultGeoManager& default_geo_manager) {
    int return_status = 0;
    this->SetIndex();
    this->SetupConnectionParameters(this->geometry_cell_type_);
    this->k0DefaultGeoShapeType_ = shape_type;
    return_status = this->GeometryInfo3DInterface::InitialGeometry(
        dx, shape_type, default_geo_manager);
    return return_status;
}
int GeometryInfoConnection3D::UpdateGeometry(
    const DefaultGeoManager& default_geo_manager) {
    return 0;
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
* @param[out]  ptr_coordi_min   minimum coordinates of the geometry.
* @param[out]  ptr_coordi_max   maximum coordinates of the geometry.
*/
void GeometryInfoConnection3D::InitialCoordinateGivenLevel(
    std::vector<DefReal>* const ptr_coordi_min,
    std::vector<DefReal>* const ptr_coordi_max) {
    vertex_given_level_.push_back({});
    connection_vertex_given_level_.push_back({});
    DefSizet i_vertex = 0;
    GeometryConnectionCoordinate vertex_temp;
    vertex_temp.coordinates = { 0., 0., 0. };
    ptr_coordi_min->resize(3);
    ptr_coordi_max->resize(3);
    ptr_coordi_min->at(kXIndex) =
        coordinate_origin_.at(0).coordinate.at(kXIndex);
    ptr_coordi_min->at(kYIndex) =
        coordinate_origin_.at(0).coordinate.at(kYIndex);
    ptr_coordi_min->at(kZIndex) =
        coordinate_origin_.at(0).coordinate.at(kZIndex);
    ptr_coordi_max->at(kXIndex) =
        coordinate_origin_.at(0).coordinate.at(kXIndex);
    ptr_coordi_max->at(kYIndex) =
        coordinate_origin_.at(0).coordinate.at(kYIndex);
    ptr_coordi_max->at(kZIndex) =
        coordinate_origin_.at(0).coordinate.at(kZIndex);
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
        if (ptr_coordi_min->at(kZIndex) >
            iter_vertex.coordinate.at(kZIndex)) {
            ptr_coordi_min->at(kZIndex) =
                iter_vertex.coordinate.at(kZIndex);
        } else if (ptr_coordi_max->at(kZIndex) <
            iter_vertex.coordinate.at(kZIndex)) {
            ptr_coordi_max->at(kZIndex) =
                iter_vertex.coordinate.at(kZIndex);
        }
        vertex_temp.coordinates.at(0) = iter_vertex.coordinate.at(0);
        vertex_temp.coordinates.at(1) = iter_vertex.coordinate.at(1);
        vertex_temp.coordinates.at(2) = iter_vertex.coordinate.at(2);
        vertex_given_level_.at(0).vec_vertex_coordinate
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
    const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
    const std::pair<DefAmrIndexUint, DefSizet>& vertex1) {
    DefReal x_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kXIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate.at(vertex1.second).coordinates.at(kXIndex));
    DefReal y_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kYIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate.at(vertex1.second).coordinates.at(kYIndex));
    DefReal z_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kZIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate.at(vertex1.second).coordinates.at(kZIndex));
    return sqrt(x_dis * x_dis + y_dis * y_dis + z_dis * z_dis);
}
/**
* @brief function to compute coordinates of the mid point of two vertices (3D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @param[out]  ptr_coordinates mid point coordinates.
*/
void GeometryInfoConnection3D::ComputeMidCoordinates(
    const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
    const std::pair<DefAmrIndexUint, DefSizet>& vertex1,
    std::vector<DefReal>* const ptr_coordinates) {
    ptr_coordinates->at(kXIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kXIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate
        .at(vertex1.second).coordinates.at(kXIndex)) / 2.;
    ptr_coordinates->at(kYIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kYIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate
        .at(vertex1.second).coordinates.at(kYIndex)) / 2.;
    ptr_coordinates->at(kZIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate.at(vertex0.second).coordinates.at(kZIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate
        .at(vertex1.second).coordinates.at(kZIndex)) / 2.;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_3D_FUNCTIONS