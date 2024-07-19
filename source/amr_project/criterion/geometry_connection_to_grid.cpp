//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
* @param[in]  i_geo   the ith geometry.
* @param[in]  i_level   level of geometry.
* @param[in]  grid_extend_type type of extension.
* @param[in]  sfbitset_aux   class manage space filling curves.
* @param[out]  ptr_grid_info grid with updated tracking node information.
*/
void GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(DefAmrIndexUint i_geo, DefAmrIndexUint i_level,
    const EGridExtendType grid_extend_type, const SFBitsetAuxInterface& sfbitset_aux,
    GridInfoInterface* const ptr_grid_info) {
    DefAmrIndexUint dims = DefAmrIndexUint(ptr_grid_info->grid_space_.size());
    if (ptr_grid_info->grid_space_.size() != vertex_given_level_.at(0)
        .vec_vertex_coordinate.at(0).coordinates.size()) {
        LogManager::LogError("Size of grid_space ("
            + std::to_string(ptr_grid_info->grid_space_.size()) +
            ") is not equal to size of coordinates ("
            + std::to_string(vertex_given_level_.at(0).vec_vertex_coordinate
            .at(0).coordinates.size()) + ") in vertex_given_level_ "
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefAmrIndexUint level_diff = ptr_grid_info->i_level_ - i_level;

    std::pair<ECriterionType, DefAmrIndexUint> key_tracking_grid = { ECriterionType::kGeometry, i_geo };
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ = grid_extend_type;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
        ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset sfbitset_tmp;
    for (auto& iter : connection_vertex_given_level_.at(level_diff)) {
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(
            ptr_grid_info->grid_space_, vertex_given_level_.at(iter.first)
            .vec_vertex_coordinate.at(iter.second).coordinates);
        vertex_given_level_.at(iter.first).vec_vertex_coordinate
            .at(iter.second).map_bitset_ref
            .insert({ ptr_grid_info->i_level_, sfbitset_tmp });
        if (ptr_grid_info->CheckIfNodeOutsideCubicDomain(dims, sfbitset_tmp, sfbitset_aux)
            >= GridInfoInterface::kFlagInsideDomain_) {
            if (ptr_tracking_node->find(sfbitset_tmp) == ptr_tracking_node->end()) {
                ptr_tracking_node->insert({ sfbitset_tmp, (ptr_grid_info
                ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())
                ->k0TrackNodeInstance_ });
            }
        } else {
            std::vector<DefReal>& coordi = vertex_given_level_.at(iter.first)
                .vec_vertex_coordinate.at(iter.second).coordinates;
            if (dims == 2) {
                LogManager::LogError("coordinate (" + std::to_string(coordi[kXIndex]) + ", "
                + std::to_string(coordi[kYIndex]) + ") is outside the computational domain in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            } else {
                LogManager::LogError("coordinate (" + std::to_string(coordi[kXIndex]) + ", "
                + std::to_string(coordi[kYIndex]) + ", "
                + std::to_string(coordi[kZIndex]) + ") is outside the computational domain in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
        }
        ptr_tracking_node->at(sfbitset_tmp).set_point_index.insert(iter);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
