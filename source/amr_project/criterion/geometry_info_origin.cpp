//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.cpp
* @author Zhengliang Liu
* @brief functions for geometries with connection relations, e.g. stl files
* @date  2023-4-25
*/
#include <vector>
#include <string>
#include "./auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_origin.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject {
void GeometryInfoOrigin::InitialGeometry(const DefReal dx) {
    this->GeometryInfoInterface::InitialGeometry(dx);
}
/**
* @brief   function to find tracking nodes based on geometries.
* @param[in] sfbitset_aux class to manage functions for space filling code computation.
* @param[in] ptr_grid_info pointer to class store grid information.
*/
void GeometryInfoOrigin::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) {
    std::pair<ECriterionType, DefInt> key_tracking_grid = { ECriterionType::kGeometry, i_geo_ };
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ = grid_extend_type_;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
     ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset sfbitset_tmp;
    std::vector<DefReal> coordi(k0GeoDim_, 0.);
    DefSizet ipoint = 0;
    for (const auto& iter : vec_vertices_) {
        coordi[kXIndex] = iter->coordinate_[kXIndex];
        coordi[kYIndex] = iter->coordinate_[kYIndex];
        if (k0GeoDim_ == 3) {
            coordi[kZIndex] = iter->coordinate_[kZIndex];
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(ptr_grid_info->GetGridSpace(), coordi);
        if (ptr_grid_info->CheckIfNodeOutsideCubicDomain(k0GeoDim_, sfbitset_tmp, sfbitset_aux)
            >= GridInfoInterface::kFlagInsideDomain_) {
            if (ptr_tracking_node->find(sfbitset_tmp) == ptr_tracking_node->end()) {
                ptr_tracking_node->insert({ sfbitset_tmp, (ptr_grid_info
                ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())->k0TrackNodeInstance_ });
            }
        } else {
            LogManager::LogError("coordinate (" + std::to_string(iter->coordinate_[kXIndex]) + ", "
                + std::to_string(iter->coordinate_[kYIndex]) + ", "
                + std::to_string(iter->coordinate_[kZIndex]) + ") is outside the computational domain");
        }
        ptr_tracking_node->at(sfbitset_tmp).set_point_index.insert({0, ipoint});
        ++ipoint;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
