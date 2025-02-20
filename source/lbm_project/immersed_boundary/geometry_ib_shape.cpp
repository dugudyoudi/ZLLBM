//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_ib_shape_2d.cpp
* @author Zhengliang Liu
* @brief functions to update 2d geometry shapes for immersed boundary method
* @date  2022-8-5
*/
#include "io/log_write.h"
#include "immersed_boundary/geometry_ib_shape.h"
namespace rootproject {
namespace lbmproject {
std::unique_ptr<amrproject::GeoShapeInterface> GeoIBShapeReader::ReadAndCreateGeoShape(
    const DefInt dims, const std::string& shape_id,
    const std::weak_ptr<amrproject::GeometryInfoInterface>& ptr_geo) const {
    if (dims == 2) {
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
        if (shape_id == "cylinder") {
            return std::make_unique<GeoShapeIBCircle2D>(ptr_geo);
        } else if (shape_id == "line") {
            return std::make_unique<GeoShapeIBLine2D>(ptr_geo);
        } else {
            amrproject::LogManager::LogError("Unsupported shape: " + shape_id);
            return nullptr;
        }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else if (dims == 3) {
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
        if (shape_id == "sphere") {
            return std::make_unique<GeoShapeIBSphere3D>(ptr_geo);
        } else if (shape_id == "cubic") {
            return std::make_unique<GeoShapeIBCubic3D>(ptr_geo);
        } else if (shape_id == "quadrilateral") {
            return std::make_unique<GeoShapeIBQuadrilateral3D>(ptr_geo);
        } else {
            amrproject::LogManager::LogError("Unsupported shape: " + shape_id);
            return nullptr;
        }
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    } else {
        amrproject::LogManager::LogError("Dimension should be 2 or 3");
        return nullptr;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
