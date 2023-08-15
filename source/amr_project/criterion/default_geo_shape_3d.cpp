//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @file default_geo_shape_3d.cpp
* @author Zhengliang Liu
* @brief functions to generate default geomtries
* @date  2023-2-13
*/
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to generate a 3D cube
*/
void DefaultGeoManager::cube_initial(
    const DefReal dx,
    GeometryInfo3DInterface* const ptr_geo) const {
    DefReal length = 0.5;
    DefAmrUint num_point = static_cast<DefAmrUint>(length / dx + kEps);
    ptr_geo->flood_fill_origin_ = ptr_geo->geometry_center_;
    ptr_geo->coordinate_origin_ = std::vector<GeometryCoordinate3D>(
        6 * num_point * num_point);
    DefAmrUint num_sum = 0;
    DefReal x_coordi, y_coordi, z_coordi;
    z_coordi = ptr_geo->geometry_center_[kZIndex] - length / 2;
    for (DefAmrUint iy = 0; iy < num_point; ++iy) {
        for (DefAmrUint ix = 0; ix < num_point; ++ix) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = dx / 2 + ix * dx - length / 2
                + ptr_geo->geometry_center_[kXIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = dx / 2 + iy * dx - length / 2
                + ptr_geo->geometry_center_[kYIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    z_coordi = ptr_geo->geometry_center_[kZIndex] + length / 2;
    for (DefAmrUint iy = 0; iy < num_point; ++iy) {
        for (DefAmrUint ix = 0; ix < num_point; ++ix) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = dx / 2 + ix * dx - length / 2
                + ptr_geo->geometry_center_[kXIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = dx / 2 + iy * dx - length / 2
                + ptr_geo->geometry_center_[kYIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = z_coordi;
            num_sum++;
        }
    }
    y_coordi = ptr_geo->geometry_center_[kYIndex] - length / 2;
    for (DefAmrUint iz = 0; iz < num_point; ++iz) {
        for (DefAmrUint ix = 0; ix < num_point; ++ix) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = dx / 2 + ix * dx - length / 2
                + ptr_geo->geometry_center_[kXIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = y_coordi;
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = dx / 2 + iz * dx - length / 2
                + ptr_geo->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    y_coordi = ptr_geo->geometry_center_[kYIndex] + length / 2;
    for (DefAmrUint iz = 0; iz < num_point; ++iz) {
        for (DefAmrUint ix = 0; ix < num_point; ++ix) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = dx / 2 + ix * dx - length / 2
                + ptr_geo->geometry_center_[kXIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = y_coordi;
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = dx / 2 + iz * dx - length / 2
                + ptr_geo->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    x_coordi = ptr_geo->geometry_center_[kXIndex] - length / 2;
    for (DefAmrUint iz = 0; iz < num_point; ++iz) {
        for (DefAmrUint iy = 0; iy < num_point; ++iy) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = x_coordi;
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = dx / 2 + iy * dx - length / 2
                + ptr_geo->geometry_center_[kYIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = dx / 2 + iz * dx - length / 2
                + ptr_geo->geometry_center_[kZIndex];
            num_sum++;
        }
    }
    x_coordi = ptr_geo->geometry_center_[kXIndex] + length / 2;
    for (DefAmrUint iz = 0; iz < num_point; ++iz) {
        for (DefAmrUint iy = 0; iy < num_point; ++iy) {
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kXIndex) = x_coordi;
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kYIndex) = dx / 2 + iy * dx - length / 2
                + ptr_geo->geometry_center_[kYIndex];
            ptr_geo->coordinate_origin_.at(num_sum)
                .coordinate.at(kZIndex) = dx / 2 + iz * dx - length / 2
                + ptr_geo->geometry_center_[kZIndex];
            num_sum++;
        }
    }
}
void DefaultGeoManager::cube_update(const DefReal sum_t,
    GeometryInfo3DInterface* const ptr_geo) const {
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
