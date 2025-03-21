//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape_3d.cpp
* @author Zhengliang Liu
* @brief functions to initialize and update 3d geometry shapes
* @date  2022-8-5
*/
#include <limits>
#include <functional>
#include <algorithm>
#include "./auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to generate a 3d cube shape.
* @param[in] dx reference spatial step.
*/
void GeoShapeDefaultCube3D::InitialShape(const DefReal dx) {
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
* @brief   function to read and set parameters for a 3d cube.
* @param[in, out] ptr_shape_parameters pointer to map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultCube3D::ReadAndSetGeoShapeParameters(
    std::map<std::string, ParserData>* const ptr_shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("cube.length", ptr_shape_parameters, &length_)) {
        LogManager::LogWarning("length is not found for 3d cubic");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("cube.center", ptr_shape_parameters, &center_)) {
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
    for (DefInt i = 0; i < edge1_num_points_; ++i) {
        mid_x = neighbor_point_.at(kXIndex) + (i + 0.5) * edge1_dir_x * edge1_arc + offset[kXIndex];
        mid_y = neighbor_point_.at(kYIndex) + (i + 0.5) * edge1_dir_y * edge1_arc + offset[kYIndex];
        mid_z = neighbor_point_.at(kZIndex) + (i + 0.5) * edge1_dir_z * edge1_arc + offset[kZIndex];
        for (DefInt j = 0; j < edge2_num_points_; ++j) {
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
* @param[in, out] ptr_shape_parameters pointer to map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultQuadrilateral3D::ReadAndSetGeoShapeParameters(
    std::map<std::string, ParserData>* const ptr_shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("quadrilateral.start_point", ptr_shape_parameters, &start_point_)) {
        LogManager::LogWarning("start_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("quadrilateral.neighbor_point", ptr_shape_parameters, &neighbor_point_)) {
        LogManager::LogWarning("neighbor_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("quadrilateral.diagonal_point", ptr_shape_parameters, &diagonal_point_)) {
        LogManager::LogWarning("diagonal_point is not found for 3d quadrilateral");
        parameter_exist = false;
    }
    return parameter_exist;
}
// simplex method to find the minimum integer, used in sphere generation with unit vectors
DefInt SimplexMethod(std::function<DefInt(DefInt)> func, DefInt n_phi_min, DefInt n_phi_max, DefInt max_iter = 100) {
    std::vector<DefInt> vertices = {n_phi_min, (n_phi_min + n_phi_max)/2, n_phi_max};
    DefReal alpha = 1.0;
    DefReal gamma = 2.0;
    DefReal rho = 0.5;
    DefReal sigma = 0.5;

    for (DefInt iter = 0; iter < max_iter; ++iter) {
        std::sort(vertices.begin(), vertices.end(), [&](DefInt a, DefInt b) {
            return func(a) < func(b);
        });

        DefInt best = vertices[0];
        DefInt worst = vertices[2];
        DefInt next_worst = vertices[1];

        DefInt centroid = (vertices[0] + vertices[1]) / 2;

        // Reflection
        DefInt reflected = centroid + static_cast<DefInt>(alpha * (centroid - worst));
        if (reflected < n_phi_min) reflected = n_phi_min;
        if (reflected > n_phi_max) reflected = n_phi_max;

        if (func(reflected) < func(next_worst)) {
            worst = reflected;
        } else {
            if (func(reflected) < func(best)) {
                DefInt expanded = centroid + static_cast<DefInt>(gamma * (reflected - centroid));
                if (expanded < n_phi_min) expanded = n_phi_min;
                if (expanded > n_phi_max) expanded = n_phi_max;

                if (func(expanded) < func(reflected)) {
                    worst = expanded;
                } else {
                    worst = reflected;
                }
            } else {
                DefInt contracted = centroid + static_cast<DefInt>(rho * (worst - centroid));
                if (contracted < n_phi_min) contracted = n_phi_min;
                if (contracted > n_phi_max) contracted = n_phi_max;

                if (func(contracted) < func(worst)) {
                    worst = contracted;
                } else {
                    for (size_t i = 1; i < vertices.size(); ++i) {
                        vertices[i] = best + static_cast<DefInt>(sigma * (vertices[i] - best));
                        if (vertices[i] < n_phi_min) vertices[i] = n_phi_min;
                        if (vertices[i] > n_phi_max) vertices[i] = n_phi_max;
                    }
                }
            }
        }
        vertices[2] = worst;

        if (vertices[0] == vertices[1] && vertices[1] == vertices[2]) {
            break;
        }
    }

    return vertices[0];
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
    if (method_generate_sphere_ == 1) {
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
    }


    //  sphere represented by unit vectors
    if (method_generate_sphere_ == 2) {
        DefInt nmax = 4;
        const DefReal epsilon = dx/radius_/2;
        const DefInt n_phi_min = static_cast<DefInt>(kPi /2/ epsilon);

        auto func_compute_total = [epsilon](DefInt n_phi) -> DefInt {
            DefInt n_total = 0, n_theta;
            DefReal phi, aa, a1, a2;
            for (DefInt i = 2; i < n_phi + 1; ++i) {
                phi = static_cast<DefReal>(i - 1) * kPi / static_cast<DefReal>(n_phi - 1) / 2;
                aa = phi + kPi / 2 / (n_phi - 1);
                a1 = cos(epsilon) - cos(phi) * cos(aa);
                a2 = sin(phi) * sin(aa);
                n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
                n_total += n_theta;
            }
            if (n_total < 1) {
                return (std::numeric_limits<DefInt>::max)();
            } else {
                return n_total;
            }
        };
        DefInt  n_phi_out = SimplexMethod(func_compute_total, n_phi_min, n_phi_min * nmax),
            min_val = func_compute_total(n_phi_out), n_theta;
        DefReal phi, aa, a1, a2;

        phi = kPi / 2;
        aa = phi + kPi / 2 / (n_phi_out - 1);
        a1 = cos(epsilon) - cos(phi) * cos(aa);
        a2 = sin(phi) * sin(aa);
        DefInt n_theta_mid = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);

        DefReal theta;
        DefInt icount = 0;
        DefInt num_points = 2*(min_val + 1) - n_theta_mid;
        ptr_geo_info->vec_vertices_.resize(num_points);
        DefInt index_tmp;
        for (DefInt j = 1; j < n_phi_out + 1; ++j) {
            phi = static_cast<DefReal>(j - 1) * kPi / static_cast<DefReal>(n_phi_out - 1) / 2;
            aa = phi + kPi / 2 / (n_phi_out - 1);
            a1 = cos(epsilon) - cos(phi) * cos(aa);
            a2 = sin(phi) * sin(aa);
            if (j == 1) {
                ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
                ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
                    + center_.at(kXIndex) + offset[kXIndex];
                ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = 0.
                    + center_.at(kYIndex) + offset[kYIndex];
                ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * radius_
                    + center_.at(kZIndex) + offset[kZIndex];
                index_tmp = num_points - 1;
                ptr_geo_info->vec_vertices_.at(index_tmp) = ptr_geo_info->GeoIndexVertexCreator();
                ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kXIndex] = -cos(phi) * radius_
                    + center_.at(kXIndex) + offset[kXIndex];
                ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kYIndex] = 0.
                    + center_.at(kYIndex) + offset[kYIndex];
                ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kZIndex] = sin(phi) * radius_
                    + center_.at(kZIndex) + offset[kZIndex];
                ++icount;
            } else if (j == n_phi_out) {
                n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
                for (DefInt k = 1; k < n_theta + 1; ++k) {
                    theta = static_cast<DefInt>(k - 1) * 2. * kPi / n_theta;
                    ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
                        + center_.at(kXIndex) + offset[kXIndex];
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
                        + center_.at(kYIndex) + offset[kYIndex];
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
                        + center_.at(kZIndex) + offset[kZIndex];
                    ++icount;
                }
            } else {
                n_theta = static_cast<DefInt>(kPi / acos(a1 / a2) + 0.5);
                for (DefInt k = 1; k < n_theta + 1; ++k) {
                    theta = static_cast<DefReal>(k - 1) * 2. * kPi / n_theta;
                    ptr_geo_info->vec_vertices_.at(icount) = ptr_geo_info->GeoIndexVertexCreator();
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kXIndex] = cos(phi) * radius_
                        + center_.at(kXIndex) + offset[kXIndex];
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
                        + center_.at(kYIndex) + offset[kYIndex];
                    ptr_geo_info->vec_vertices_.at(icount)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
                        + center_.at(kZIndex) + offset[kZIndex];
                    index_tmp = num_points - icount - 1;
                    ptr_geo_info->vec_vertices_.at(index_tmp) = ptr_geo_info->GeoIndexVertexCreator();
                    ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kXIndex] = -cos(phi) * radius_
                        + center_.at(kXIndex) + offset[kXIndex];
                    ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kYIndex] = sin(phi) * sin(theta) * radius_
                        + center_.at(kYIndex) + offset[kYIndex];
                    ptr_geo_info->vec_vertices_.at(index_tmp)->coordinate_[kZIndex] = sin(phi) * cos(theta) * radius_
                        + center_.at(kZIndex) + offset[kZIndex];
                    ++icount;
                }
            }
        }
    }
}
/**
* @brief   function to read and set parameters for a 3d sphere.
* @param[in, out] ptr_shape_parameters pointer to map storing shape parameters.
* @return true if all parameters are found, false if some parameters are missing.
*/
bool GeoShapeDefaultSphere3D::ReadAndSetGeoShapeParameters(
    std::map<std::string, ParserData>* const ptr_shape_parameters) {
    InputParser input_parser;
    bool parameter_exist = true;
    if (!input_parser.GetValue<DefReal>("sphere.center", ptr_shape_parameters, &center_)) {
        LogManager::LogWarning("center is not found for 3d sphere");
        parameter_exist = false;
    }
    if (!input_parser.GetValue<DefReal>("sphere.radius", ptr_shape_parameters, &radius_)) {
        LogManager::LogWarning("neighbor_point is not found for 3d sphere");
        parameter_exist = false;
    }
    return parameter_exist;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
