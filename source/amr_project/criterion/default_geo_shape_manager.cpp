//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file default_geo_shape_manager.cpp
* @author Zhengliang Liu
* @brief functions to generate default geomtries
* @date  2022-8-5
*/
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
namespace criterion {
/**
* @brief   function to generate a 2D circle
*/
void DefaultGeoShapeManager::geo_circle_initial(const DefUint dims,
    std::shared_ptr<GeometryInfo2DInterface> geo_ptr) {
    geo_ptr->geometry_cell_type_ = EGeometryCellType::kPolyLine;
    DefReal radius = 0.5;
    DefSizet num_vertexs = 100;
    geo_ptr->flood_fill_origin_ = geo_ptr->geometry_center_;
    geo_ptr->coordinate_origin_ = std::vector<GeometryCoordinate2D>(num_vertexs);
    if (dims == 2) {
        for (DefSizet i = 0; i < num_vertexs; ++i) {
            geo_ptr->coordinate_origin_.at(i).coordinate[kXIndex] = radius
                * cos(2. * kPi * (static_cast<DefReal>(i)
                    / static_cast<DefReal>(num_vertexs)))
                + geo_ptr->geometry_center_[kXIndex]
                + geo_ptr->k0RealOffest_[kXIndex];
            geo_ptr->coordinate_origin_.at(i).coordinate[kYIndex] = radius
                * sin(2. * kPi * (static_cast<DefReal>(i)
                    / static_cast<DefReal>(num_vertexs)))
                + geo_ptr->geometry_center_[kYIndex]
                + geo_ptr->k0RealOffest_[kYIndex];
        }
    }
}
void DefaultGeoShapeManager::geo_circle_update(const DefUint dims,
    std::shared_ptr<GeometryInfo2DInterface> geo_ptr) {
}
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
