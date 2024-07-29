//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.cpp
* @author Zhengliang Liu
* @brief functions to manage geometry information
* @date  2022-8-25
*/
#include "io/log_write.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to initialize geometry.
 * @param dx reference spatial step.
 */
int GeometryInfoInterface::InitialGeometry(const DefReal dx) {
    ptr_geo_shape_->ptr_geo_info_ = this;
    ptr_geo_shape_->InitialShape(dx);
    return 0;
}
/**
 * @brief function to update geometry.
 * @param sum_t time step at current level.
 * @return 0 the geometry is updated, 1 the geometry is not updated
 */
int GeometryInfoInterface::UpdateGeometry(const DefReal sum_t) {
    if (geometry_status_ == EGeometryStatus::kMoving) {
        ptr_geo_shape_->UpdateShape(sum_t);
    } else {
        return -1;
    }
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
