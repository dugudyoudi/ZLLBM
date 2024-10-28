//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary_geometry.h
* @author Zhengliang Liu
* @brief define classes to manage geometries used in immersed boundary method.
* @date  2024-10-17
*/
#ifndef SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_GEOMETRY_H_
#define SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_GEOMETRY_H_
#include <vector>
#include <array>
#include <map>
#include <memory>
#include "criterion/geometry_info_origin.h"
#include "immersed_boundary/immersed_boundary.h"
namespace rootproject {
namespace lbmproject {
class GeometryInfoImmersedBoundary : public amrproject::GeometryInfoOrigin, public FsiImmersedBoundary {
 public:
    std::unique_ptr<amrproject::GeometryVertex> GeoVertexCreator() const override {
        return std::make_unique<GeometryVertexImmersedBoundary>();
    }
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_GEOMETRY_H_
