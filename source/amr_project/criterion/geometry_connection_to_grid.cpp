//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_connection_to_grid.cpp
* @author Zhengliang Liu
* @brief functions for grid generation based on geometries with connection relations
* @date  2022-10-20
*/
#include <algorithm>
#include <queue>
#include <vector>
#include "./auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_connection.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to find tracking near the geometry described by connection relations.
* @param[in]  dims   dimension of the geometry.
* @param[in]  i_geo   the ith geometry.
* @param[in]  i_level   level of geometry.
* @param[in]  grid_extend_type type of extension.
* @param[in]  sfbitset_aux   class manage space filling curves.
* @param[out]  ptr_grid_info grid with updated tracking node information.
*/
void GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(DefInt dims,
    DefInt i_geo, DefInt i_level,
    const EGridExtendType grid_extend_type, const SFBitsetAuxInterface& sfbitset_aux,
    GridInfoInterface* const ptr_grid_info) {
    const DefInt i_level_grid = ptr_grid_info->GetGridLevel();
    const DefInt level_diff = i_level_grid - i_level;
    const std::vector<DefReal>& grid_space = ptr_grid_info->GetGridSpace();

    std::pair<ECriterionType, DefInt> key_tracking_grid = { ECriterionType::kGeometry, i_geo };
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ = grid_extend_type;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
        ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset sfbitset_tmp;
    for (auto& iter : connection_vertex_given_level_.at(level_diff)) {
        std::array<DefReal, 3>& coordinate = vertex_given_level_.at(iter.first)
            .vec_vertex_coordinate_.at(iter.second)->coordinate_;
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(grid_space,
            {coordinate[kXIndex], coordinate[kYIndex], coordinate[kZIndex]});
        vertex_given_level_.at(iter.first).vec_vertex_coordinate_
            .at(iter.second)->map_bitset_ref_.insert({i_level_grid, sfbitset_tmp});
        if (ptr_grid_info->CheckIfNodeOutsideCubicDomain(dims, sfbitset_tmp, sfbitset_aux)
            >= GridInfoInterface::kFlagInsideDomain_) {
            if (ptr_tracking_node->find(sfbitset_tmp) == ptr_tracking_node->end()) {
                ptr_tracking_node->insert({ sfbitset_tmp, (ptr_grid_info
                ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())
                ->k0TrackNodeInstance_ });
            }
        } else {
            std::array<DefReal, 3>& coordi = vertex_given_level_.at(iter.first)
                .vec_vertex_coordinate_.at(iter.second)->coordinate_;
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
