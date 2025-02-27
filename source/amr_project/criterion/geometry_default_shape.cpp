//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape.cpp
* @author Zhengliang Liu
* @brief functions to manage default geometry shapes
*/
#include <limits>
#include <memory>
#include "io/log_write.h"
#include "io/input_parser.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
std::unique_ptr<GeoShapeInterface> GeoShapeReader::ReadAndCreateGeoShape(const DefInt dims,
    const std::string& shape_id, const std::weak_ptr<GeometryInfoInterface>& ptr_geo) const {

    if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
        if (shape_id == "cylinder") {
            return std::make_unique<GeoShapeDefaultCircle2D>(ptr_geo);
        } else if (shape_id == "line") {
            return std::make_unique<GeoShapeDefaultLine2D>(ptr_geo);
        } else {
            LogManager::LogError("Unsupported shape: " + shape_id);
            return nullptr;
        }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else if (dims == 3) {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
        if (shape_id == "sphere") {
            return std::make_unique<GeoShapeDefaultSphere3D>(ptr_geo);
        } else if (shape_id == "cube") {
            return std::make_unique<GeoShapeDefaultCube3D>(ptr_geo);
        } else if (shape_id == "quadrilateral") {
            return std::make_unique<GeoShapeDefaultQuadrilateral3D>(ptr_geo);
        } else {
            LogManager::LogError("Unsupported shape: " + shape_id);
            return nullptr;
        }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else {
        LogManager::LogError("Dimension should be 2 or 3");
        return nullptr;
    }
}

}  // end namespace amrproject
}  // end namespace rootproject
