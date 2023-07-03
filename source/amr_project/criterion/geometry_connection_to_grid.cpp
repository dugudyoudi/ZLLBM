//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.cpp
* @author Zhengliang Liu
* @brief functions for geometries with connection relations, e.g. stl files
* @date  2022-10-20
*/
#include <algorithm>
#include <queue>
#include <vector>
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_connection.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to find tracking near the geometry described
*               by connection relations.
* @param[in]  ptr_sfbitset_aux   pointer to class manage space filling curves.
* @param[out]  ptr_grid_info grid with updated tracking node information.
*/
void GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface* ptr_sfbitset_aux,
    GridInfoInterface* const ptr_grid_info) {
    if (ptr_grid_info->grid_space_.size() != vertex_given_level_.at(0)
        .vec_vertex_coordinate.at(0).coordinates.size()) {
        LogError("Size of grid_space ("
            + std::to_string(ptr_grid_info->grid_space_.size()) +
            ") is not equal to size of coordinates("
        + std::to_string(vertex_given_level_.at(0).vec_vertex_coordinate
            .at(0).coordinates.size()) + ") in vertex_given_level_.");
    }
    DefSizet level_diff = ptr_grid_info->i_level_ - i_level_;

    // create instance of tracking grid for the given geometry
    std::pair<ECriterionType, DefSizet> key_tracking_grid =
    { ECriterionType::kGeometry, i_geo_ };
    if (ptr_grid_info->map_ptr_tracking_grid_info_.find(key_tracking_grid)
        == ptr_grid_info->map_ptr_tracking_grid_info_.end()) {
        ptr_grid_info->map_ptr_tracking_grid_info_.insert({ key_tracking_grid,
            ptr_tracking_grid_info_creator_->CreateTrackingGridInfo() });
    }
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ =
        grid_extend_type_;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
        ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()
        ->map_tracking_node_);

    DefSFBitset bitset_temp;
    for (auto& iter : connection_vertex_given_level_.at(level_diff)) {
        bitset_temp = ptr_sfbitset_aux->SFBitsetEncodingCoordi(
            ptr_grid_info->grid_space_, vertex_given_level_.at(iter.first)
            .vec_vertex_coordinate.at(iter.second).coordinates);
        vertex_given_level_.at(iter.first).vec_vertex_coordinate
            .at(iter.second).map_bitset_ref
            .insert({ ptr_grid_info->i_level_, bitset_temp });
        if (ptr_tracking_node->find(bitset_temp) == ptr_tracking_node->end()) {
            ptr_tracking_node->insert({ bitset_temp, (ptr_grid_info
            ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())
            ->k0TrackNodeInstance_ });
        }
        ptr_tracking_node->at(bitset_temp).set_point_index.insert(iter);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
