//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape_3d.cpp
* @author Zhengliang Liu
* @brief functions to initialize and update 3d geometry shapes
* @date  2022-8-5
*/
#include <limits>
#include "./auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to generate a 3d cubic shape.
* @param[in] dx reference spatial step.
*/
void GeoShapeDefaultCubic3D::InitialShape(const DefReal dx) {
    DefInt num_point_dir = static_cast<DefInt>(length_ / dx + kEps);
    GeometryInfoInterface* ptr_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_geo_info = ptr.get();
    } else {
        LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    ptr_geo_info->SetGeometryCenter(center_);
    const std::array<DefReal, 3>& center = ptr_geo_info->GetGeometryCenter();
    ptr_geo_info->SetFloodFillOrigin(center);
    DefSizet num_points = 6 * num_point_dir * num_point_dir;
    ptr_geo_info->vec_vertices_.resize(num_points);
    DefInt num_sum = 0;
    DefReal x_coordi, y_coordi, z_coordi;
    z_coordi = center[kZIndex] - length_ / 2;
    for (DefInt iy = 0; iy < num_point_dir; ++iy) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kXIndex) = dx / 2 + ix * dx - length_ / 2 + center[kXIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = dx / 2 + iy * dx - length_ / 2 + center[kYIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    z_coordi = center[kZIndex] + length_ / 2;
    for (DefInt iy = 0; iy < num_point_dir; ++iy) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kXIndex) = dx / 2 + ix * dx - length_ / 2 + center[kXIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = dx / 2 + iy * dx - length_ / 2 + center[kYIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    y_coordi = center[kYIndex] - length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kXIndex) = dx / 2 + ix * dx - length_ / 2 + center[kXIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = y_coordi;
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = dx / 2 + iz * dx - length_ / 2 + center[kZIndex];
            num_sum++;
        }
    }
    y_coordi = center[kYIndex] + length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kXIndex) = dx / 2 + ix * dx - length_ / 2 + center[kXIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = y_coordi;
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = dx / 2 + iz * dx - length_ / 2 + center[kZIndex];
            num_sum++;
        }
    }
    x_coordi = center[kXIndex] - length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt iy = 0; iy < num_point_dir; ++iy) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)->coordinate_.at(kXIndex) = x_coordi;
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = dx / 2 + iy * dx - length_ / 2 + center[kYIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = dx / 2 + iz * dx - length_ / 2 + center[kZIndex];
            num_sum++;
        }
    }
    x_coordi = center[kXIndex] + length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt iy = 0; iy < num_point_dir; ++iy) {
            ptr_geo_info->vec_vertices_.at(num_sum) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(num_sum)->coordinate_.at(kXIndex) = x_coordi;
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kYIndex) = dx / 2 + iy * dx - length_ / 2 + center[kYIndex];
            ptr_geo_info->vec_vertices_.at(num_sum)
                ->coordinate_.at(kZIndex) = dx / 2 + iz * dx - length_ / 2 + center[kZIndex];
            num_sum++;
        }
    }
}
/**
* @brief   function to read and set parameters for a 3d cubic.
* @param[in] shape_parameters map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultCubic3D::ReadAndSetGeoShapeParameters(const std::map<std::string, std::string>& shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("cubic.length", shape_parameters, &length_)) {
        LogManager::LogWarning("length is not found for 3d cubic");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("cubic.center", shape_parameters, &center_)) {
        LogManager::LogWarning("center is not found for 3d cubic");
        parameter_exist = false;
    }
    return parameter_exist;
}
/**
* @brief   function to generate a 3d quadrilateral.
* @param[in] dx_background reference spatial step (background).
*/
// (start)
// *                *
// *                *
// (neighbor)       (diagonal)
void GeoShapeDefaultQuadrilateral3D::InitialShape(const DefReal dx_background) {
    GeometryInfoInterface* ptr_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_geo_info = ptr.get();
    } else {
        LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const DefReal dx = dx_background/TwoPowerN(ptr_geo_info->GetLevel());
    std::array<DefReal, 3> offset = ptr_geo_info->GetOffset();
    ptr_geo_info->SetGeometryCenter({
        0.5*(start_point_.at(kXIndex) + diagonal_point_.at(kXIndex)) + offset[kXIndex],
        0.5*(start_point_.at(kYIndex) + diagonal_point_.at(kYIndex)) + offset[kYIndex],
        0.5*(start_point_.at(kZIndex) + diagonal_point_.at(kZIndex)) + offset[kZIndex]});
    ptr_geo_info->SetFloodFillOrigin(ptr_geo_info->GetGeometryCenter());
    DefReal edge1_length = std::sqrt(Square(start_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(start_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    DefReal edge2_length = std::sqrt(Square(diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    edge1_num_points_ = DefInt(edge1_length / dx + kEps) + 1;
    edge2_num_points_ = DefInt(edge2_length / dx + kEps) + 1;
    DefInt num_points = edge1_num_points_ * edge2_num_points_;
    ptr_geo_info->vec_vertices_.resize(num_points);
    DefReal edge1_arc = edge1_length / edge1_num_points_,
        edge2_arc = edge2_length / edge2_num_points_,
        edge1_dir_x = (start_point_.at(kXIndex) - neighbor_point_.at(kXIndex)) / edge1_length,
        edge1_dir_y = (start_point_.at(kYIndex) - neighbor_point_.at(kYIndex)) / edge1_length,
        edge1_dir_z = (start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)) / edge1_length,
        edge2_dir_x = (diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex)) / edge2_length,
        edge2_dir_y = (diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex)) / edge2_length,
        edge2_dir_z = (diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)) / edge2_length;
    DefReal mid_x, mid_y, mid_z;
    DefSizet node_index;
    for (DefSizet i = 0; i < edge1_num_points_; ++i) {
        mid_x = neighbor_point_.at(kXIndex) + (i + 0.5) * edge1_dir_x * edge1_arc + offset[kXIndex];
        mid_y = neighbor_point_.at(kYIndex) + (i + 0.5) * edge1_dir_y * edge1_arc + offset[kYIndex];
        mid_z = neighbor_point_.at(kZIndex) + (i + 0.5) * edge1_dir_z * edge1_arc + offset[kZIndex];
        for (DefSizet j = 0; j < edge2_num_points_; ++j) {
            node_index = i*edge1_num_points_ + j;
            ptr_geo_info->vec_vertices_.at(node_index) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(node_index)->coordinate_[kXIndex] =
                mid_x + (j + 0.5) * edge2_dir_x * edge2_arc;
            ptr_geo_info->vec_vertices_.at(node_index)->coordinate_[kYIndex] =
                mid_y + (j + 0.5) * edge2_dir_y * edge2_arc;
            ptr_geo_info->vec_vertices_.at(node_index)->coordinate_[kZIndex] =
                mid_z + (j + 0.5) * edge2_dir_z * edge2_arc;
        }
    }
}
/**
* @brief   function to read and set parameters for a 3d quadrilateral.
* @param[in] shape_parameters map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultQuadrilateral3D::ReadAndSetGeoShapeParameters(
    const std::map<std::string, std::string>& shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("quadrilateral.start_point", shape_parameters, &start_point_)) {
        LogManager::LogWarning("start_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("quadrilateral.neighbor_point", shape_parameters, &neighbor_point_)) {
        LogManager::LogWarning("neighbor_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("quadrilateral.diagonal_point", shape_parameters, &diagonal_point_)) {
        LogManager::LogWarning("diagonal_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    return parameter_exist;
}
/**
* @brief   function to generate a sphere
* @param[in] dx_background reference spatial step (background)
*/
void GeoShapeDefaultSphere3D::InitialShape(const DefReal dx_background) {
    GeometryInfoInterface* ptr_geo_info = nullptr;
    if (auto ptr = ptr_geo_info_.lock()) {
        ptr_geo_info = ptr.get();
    } else {
        LogManager::LogError("point to geometry information is not found for 2d circle");
    }
    const DefReal dx = dx_background/TwoPowerN(ptr_geo_info->GetLevel());
    std::array<DefReal, 3> offset = ptr_geo_info->GetOffset();
    ptr_geo_info->SetGeometryCenter({ center_.at(kXIndex) + offset[kXIndex],
        center_.at(kYIndex) + offset[kYIndex], center_.at(kZIndex) + offset[kZIndex]});
    ptr_geo_info->SetFloodFillOrigin(ptr_geo_info->GetGeometryCenter());

    // Fibonacci sphere
    const DefReal phi = kPi * (sqrt(5.) - 1.);
    const DefInt num_points = static_cast<DefInt>(4.*kPi*Square(radius_)/Square(dx) + kEps) + 1;
    ptr_geo_info->vec_vertices_.resize(num_points);
    DefReal unit_y, unit_x, unit_z, unit_r, theta;
    for (DefInt i = 0; i < num_points; ++i) {
        unit_x = 1. - (i / static_cast<DefReal>(num_points - 1)) * 2;
        unit_r = sqrt(1. - unit_x*unit_x);
        theta = phi * i;
        unit_y = cos(theta) * unit_r;
        unit_z = sin(theta) * unit_r;
        ptr_geo_info->vec_vertices_.at(i) = ptr_geo_info->GeoIndexVertexCreator();
            ptr_geo_info->vec_vertices_.at(i)->coordinate_[kXIndex] =
                unit_x*radius_ + center_.at(kXIndex) + offset[kXIndex];
            ptr_geo_info->vec_vertices_.at(i)->coordinate_[kYIndex] =
                unit_y*radius_ + center_.at(kYIndex) + offset[kYIndex];
            ptr_geo_info->vec_vertices_.at(i)->coordinate_[kZIndex] =
                unit_z*radius_ + center_.at(kZIndex) + offset[kZIndex];
    }


    // DefInt nmax = 4;
    // const DefReal epsilon = dx/(4*radius_);
    // const DefInt n_phi_min = static_cast<DefInt>(kPi / 2. / epsilon) + 2;
    // DefInt min_val = 0, n_phi_out = 0;
    // min_val = ~min_val;
    // DefReal phi, aa, a1, a2;
    // DefInt n_theta = 0, n_theta_mid = 0;
    // for (DefInt n_phi = n_phi_min; n_phi < n_phi_min * nmax + 1; ++n_phi) {
    //     DefInt  n_total = 0;
    //     for (DefInt i = 2; i < n_phi + 1; ++i) {
    //         phi = static_cast<DefReal>(i - 1) * kPi / static_cast<DefReal>(n_phi - 1) / 2;
    //         aa = phi + kPi / 2 / (n_phi - 1);
    //         a1 = cos(epsilon) - cos(phi) * cos(aa);
    //         a2 = sin(phi) * sin(aa);
    //         n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
    //         n_total += n_theta;
    //     }
    //     if (min_val > n_total) {
    //         min_val = n_total;
    //         n_phi_out = n_phi;
    //         n_theta_mid = n_theta;
    //     }
    // }
    // DefReal theta, dphi = kPi / (n_phi_out - 1) / 4.;
    // DefInt icount = 0;
    // DefInt num_points = 2*(min_val + 1) - n_theta_mid;
    // ptr_geo_info->vec_vertices_.resize(num_points);
    // DefInt index_tmp;
    // for (DefInt j = 1; j < n_phi_out + 1; ++j) {
    //     phi = static_cast<DefReal>(j - 1) * kPi / static_cast<DefReal>(n_phi_out - 1) / 2;
    //     aa = phi + kPi / 2 / (n_phi_out - 1);
    //     a1 = cos(epsilon) - cos(phi) * cos(aa);
    //     a2 = sin(phi) * sin(aa);
    //     if (j == 1) {
    //         ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
    //         ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
    //             + center_.at(kXIndex) + offset[kXIndex];
    //         ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = 0.
    //             + center_.at(kYIndex) + offset[kYIndex];
    //         ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * radius_
    //             + center_.at(kZIndex) + offset[kZIndex];
    //         index_tmp = num_points - 1;
    //         ptr_geo_info->vec_vertices_.at(index_tmp) = ptr_geo_info->GeoIndexVertexCreator();
    //         ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kXIndex] = -cos(phi) * radius_
    //             + center_.at(kXIndex) + offset[kXIndex];
    //         ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kYIndex] = 0.
    //             + center_.at(kYIndex) + offset[kYIndex];
    //         ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kZIndex] = sin(phi) * radius_
    //             + center_.at(kZIndex) + offset[kZIndex];
    //         ++icount;
    //     } else if (j == n_phi_out) {
    //         n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
    //         for (DefInt k = 1; k < n_theta + 1; ++k) {
    //             theta = static_cast<DefInt>(k - 1) * 2. * kPi / n_theta;
    //             ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
    //                 + center_.at(kXIndex) + offset[kXIndex];
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
    //                 + center_.at(kYIndex) + offset[kYIndex];
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
    //                 + center_.at(kZIndex) + offset[kZIndex];
    //             ++icount;
    //         }
    //     } else {
    //         n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
    //         for (DefInt k = 1; k < n_theta + 1; ++k) {
    //             theta = static_cast<DefReal>(k - 1) * 2. * kPi / n_theta;
    //             ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
    //                 + center_.at(kXIndex) + offset[kXIndex];
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
    //                 + center_.at(kYIndex) + offset[kYIndex];
    //             ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
    //                 + center_.at(kZIndex) + offset[kZIndex];
    //             index_tmp = num_points - icount - 1;
    //             ptr_geo_info->vec_vertices_.at(index_tmp) = ptr_geo_info->GeoIndexVertexCreator();
    //             ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kXIndex] = -cos(phi) * radius_
    //                 + center_.at(kXIndex) + offset[kXIndex];
    //             ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
    //                 + center_.at(kYIndex) + offset[kYIndex];
    //             ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
    //                 + center_.at(kZIndex) + offset[kZIndex];
    //             ++icount;
    //         }
    //     }
    // }
}
/**
* @brief   function to read and set parameters for a 3d sphere.
* @param[in] shape_parameters map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultSphere3D::ReadAndSetGeoShapeParameters(
    const std::map<std::string, std::string>& shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("sphere.center", shape_parameters, &center_)) {
        LogManager::LogWarning("center is not found for 3d sphere");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("sphere.radius", shape_parameters, &radius_)) {
        LogManager::LogWarning("neighbor_point is not found for 3d sphere");
        parameter_exist = false;
    }
    return parameter_exist;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
