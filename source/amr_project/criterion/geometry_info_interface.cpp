//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.cpp
* @author Zhengliang Liu
* @brief functions to initialize and update geometry information
* @date  2022-8-5
*/
#include <limits>
#include "io/log_write.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
int GeometryInfo2DInterface::InitialGeometry(const DefReal dx,
    const DefaultGeoShapeType shape_type, const DefaultGeoManager& default_geo_manager) {
    switch (shape_type) {
    case DefaultGeoShapeType::kCircle:
        default_geo_manager.circle_initial(this);
        return 0;
    default:
        return 1;
    }
}
int GeometryInfo2DInterface::UpdateGeometry(
    const DefaultGeoManager& default_geo_manager) {
    return 0;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
int GeometryInfo3DInterface::InitialGeometry(const DefReal dx,
    const DefaultGeoShapeType shape_type,
    const DefaultGeoManager& default_geo_manager) {
    switch (shape_type) {
    case DefaultGeoShapeType::kCube:
        default_geo_manager.cube_initial(dx, this);
        return 0;
    default:
        return 1;
    }
}
int GeometryInfo3DInterface::UpdateGeometry(
    const DefaultGeoManager& default_geo_manager) {
    return 0;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
