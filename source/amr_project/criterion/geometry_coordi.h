//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_coordi.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define struct to store geometry coordinates.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_COORDI_H_
#define ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_COORDI_H_
#include <array>
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @struct GeometryCoordinate2D
* @brief struct used to store information of a geometry vertex
*/
struct GeometryCoordinate2D {
    std::array<DefReal, 2> coordinate{};

    // GeometryCoordinate2D& operator=(const GeometryCoordinate2D& coordi_r) {
    //    this->coordinate.at(0) = coordi_r.coordinate.at(0);
    //    this->coordinate.at(1) = coordi_r.coordinate.at(1);
    //    return *this;
    // }
    // GeometryCoordinate2D operator+(const GeometryCoordinate2D& coordi_r) {
    //    GeometryCoordinate2D coordi_return;
    //    coordi_return.coordinate.at(0) = this->coordinate.at(0)
    //        + coordi_r.coordinate.at(0);
    //    coordi_return.coordinate.at(1) = this->coordinate.at(1)
    //        + coordi_r.coordinate.at(1);
    //    return coordi_return;
    // }
    // GeometryCoordinate2D operator/(DefReal real_r) {
    //    GeometryCoordinate2D coordi_return;
    //    coordi_return.coordinate.at(0) = this->coordinate.at(0) / real_r;
    //    coordi_return.coordinate.at(1) = this->coordinate.at(1) / real_r;
    // }
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @struct GeometryCoordinate3D
* @brief struct used to store information of a geometry vertex
*/
struct GeometryCoordinate3D {
    std::array<DefReal, 3> coordinate{};
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_COORDI_H_
