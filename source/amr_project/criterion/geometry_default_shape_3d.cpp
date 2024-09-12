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
* @brief   function to generate a 3d cubic shape
* @param[in] dx reference spatial step
*/
void GeoShapeDefaultCubic3D::InitialShape(const DefReal dx) {
    if (ptr_geo_info_ == nullptr) {
        LogManager::LogError("pointer to geometry infomation instance is nullptr");
    }
    DefInt num_point_dir = static_cast<DefInt>(length_ / dx + kEps);
    ptr_geo_info_->flood_fill_origin_ = ptr_geo_info_->geometry_center_;
    DefSizet num_points = 6 * num_point_dir * num_point_dir;
    ptr_geo_info_->vec_vertices_.resize(num_points);
    DefInt num_sum = 0;
    DefReal x_coordi, y_coordi, z_coordi;
    z_coordi = ptr_geo_info_->geometry_center_[kZIndex] - length_ / 2;
    for (DefInt iy = 0; iy < num_point_dir; ++iy) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kXIndex) = dx / 2 + ix * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kXIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = dx / 2 + iy * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kYIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    z_coordi = ptr_geo_info_->geometry_center_[kZIndex] + length_ / 2;
    for (DefInt iy = 0; iy < num_point_dir; ++iy) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kXIndex) = dx / 2 + ix * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kXIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = dx / 2 + iy * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kYIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    y_coordi = ptr_geo_info_->geometry_center_[kYIndex] - length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kXIndex) = dx / 2 + ix * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kXIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = y_coordi;
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = dx / 2 + iz * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    y_coordi = ptr_geo_info_->geometry_center_[kYIndex] + length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt ix = 0; ix < num_point_dir; ++ix) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kXIndex) = dx / 2 + ix * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kXIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = y_coordi;
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = dx / 2 + iz * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    x_coordi = ptr_geo_info_->geometry_center_[kXIndex] - length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt iy = 0; iy < num_point_dir; ++iy) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)->coordinate.at(kXIndex) = x_coordi;
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = dx / 2 + iy * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kYIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = dx / 2 + iz * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    x_coordi = ptr_geo_info_->geometry_center_[kXIndex] + length_ / 2;
    for (DefInt iz = 0; iz < num_point_dir; ++iz) {
        for (DefInt iy = 0; iy < num_point_dir; ++iy) {
            ptr_geo_info_->vec_vertices_.at(num_sum) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(num_sum)->coordinate.at(kXIndex) = x_coordi;
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kYIndex) = dx / 2 + iy * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kYIndex];
            ptr_geo_info_->vec_vertices_.at(num_sum)
                ->coordinate.at(kZIndex) = dx / 2 + iz * dx - length_ / 2
                + ptr_geo_info_->geometry_center_[kZIndex];
            num_sum++;
        }
    }
}
void GeoShapeDefaultCubic3D::UpdateShape(const DefReal dx) {
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to generate a quadrilateral
* @param[in] dx_background reference spatial step (background)
*/
// (start)
// *                *
// *                *
// (neighbor)       (diagonal)
void GeoShapeDefaultQuadrilateral3D::InitialShape(const DefReal dx_background) {
    if (ptr_geo_info_ == nullptr) {
        LogManager::LogError("pointer to geometry infomation instance is nullptr");
    }
    const DefReal dx = dx_background/TwoPowerN(ptr_geo_info_->i_level_);
    ptr_geo_info_->geometry_center_ = {
        start_point_.at(kXIndex) + diagonal_point_.at(kXIndex) + ptr_geo_info_->k0RealMin_[kXIndex],
        start_point_.at(kYIndex) + diagonal_point_.at(kYIndex) + ptr_geo_info_->k0RealMin_[kYIndex],
        start_point_.at(kZIndex) + diagonal_point_.at(kZIndex) + ptr_geo_info_->k0RealMin_[kZIndex]};
    ptr_geo_info_->flood_fill_origin_ = ptr_geo_info_->geometry_center_;
    DefReal edge1_length = std::sqrt(Square(start_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(start_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    DefReal edge2_length = std::sqrt(Square(diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex))
        + Square(diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex))
        + Square(diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)));
    DefSizet edge1_num_points = DefSizet(edge1_length / dx + kEps) + 1;
    DefSizet edge2_num_points = DefSizet(edge2_length / dx + kEps) + 1;
    DefSizet num_points = edge1_num_points * edge2_num_points;
    ptr_geo_info_->vec_vertices_.resize(num_points);
    DefReal edge1_arc = edge1_length / edge1_num_points,
        edge2_arc = edge2_length / edge2_num_points,
        edge1_dir_x = (start_point_.at(kXIndex) - neighbor_point_.at(kXIndex)) / edge1_length,
        edge1_dir_y = (start_point_.at(kYIndex) - neighbor_point_.at(kYIndex)) / edge1_length,
        edge1_dir_z = (start_point_.at(kZIndex) - neighbor_point_.at(kZIndex)) / edge1_length,
        edge2_dir_x = (diagonal_point_.at(kXIndex) - neighbor_point_.at(kXIndex)) / edge2_length,
        edge2_dir_y = (diagonal_point_.at(kYIndex) - neighbor_point_.at(kYIndex)) / edge2_length,
        edge2_dir_z = (diagonal_point_.at(kZIndex) - neighbor_point_.at(kZIndex)) / edge2_length;
    DefReal mid_x, mid_y, mid_z;
    DefSizet node_index;
    for (DefSizet i = 0; i < edge1_num_points; ++i) {
        mid_x = neighbor_point_.at(kXIndex) + (i + 0.5) * edge1_dir_x * edge1_arc
            + ptr_geo_info_->k0RealMin_[kXIndex];
        mid_y = neighbor_point_.at(kYIndex) + (i + 0.5) * edge1_dir_y * edge1_arc
            + ptr_geo_info_->k0RealMin_[kYIndex];
        mid_z = neighbor_point_.at(kZIndex) + (i + 0.5) * edge1_dir_z * edge1_arc
            + ptr_geo_info_->k0RealMin_[kZIndex];
        for (DefSizet j = 0; j < edge2_num_points; ++j) {
            node_index = i*edge1_num_points + j;
            ptr_geo_info_->vec_vertices_.at(node_index) = ptr_geo_info_->GeoVertexCreator();
            ptr_geo_info_->vec_vertices_.at(node_index)->coordinate[kXIndex] =
                mid_x + (j + 0.5) * edge2_dir_x * edge2_arc;
            ptr_geo_info_->vec_vertices_.at(node_index)->coordinate[kYIndex] =
                mid_y + (j + 0.5) * edge2_dir_y * edge2_arc;
            ptr_geo_info_->vec_vertices_.at(node_index)->coordinate[kZIndex] =
                mid_z + (j + 0.5) * edge2_dir_z * edge2_arc;
        }
    }
}
void GeoShapeDefaultQuadrilateral3D::UpdateShape(const DefReal sum_t) {
}
}  // end namespace amrproject
}  // end namespace rootproject
