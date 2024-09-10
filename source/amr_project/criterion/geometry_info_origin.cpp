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
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_origin.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
*/
void GeometryInfoOrigin2D::SetIndex() {
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachVertex_);
}
int GeometryInfoOrigin2D::InitialGeometry(const DefReal dx) {
    int return_status = 0;
    this->SetIndex();
    return_status = this->GeometryInfoInterface::InitialGeometry(dx);
    return return_status;
}
/**
* @brief   function to find tracking nodes based on geometries.
* @param[in] sfbitset_aux class to manage functions for space filling code computation.
* @param[in] ptr_grid_info pointer to class store grid information.
*/
void GeometryInfoOrigin2D::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) {
    std::pair<ECriterionType, DefInt> key_tracking_grid = { ECriterionType::kGeometry, i_geo_ };
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ = grid_extend_type_;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
     ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset sfbitset_tmp;
    std::vector<DefReal> coordi(2, 0.);
    DefSizet ipoint = 0;
    for (const auto& iter : coordinate_origin_) {
        coordi[kXIndex] = iter.coordinate[kXIndex];
        coordi[kYIndex] = iter.coordinate[kYIndex];
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(ptr_grid_info->grid_space_, coordi);
        if (ptr_grid_info->CheckIfNodeOutsideCubicDomain(2, sfbitset_tmp, sfbitset_aux)
            >= GridInfoInterface::kFlagInsideDomain_) {
            if (ptr_tracking_node->find(sfbitset_tmp) == ptr_tracking_node->end()) {
                ptr_tracking_node->insert({ sfbitset_tmp, (ptr_grid_info
                ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())->k0TrackNodeInstance_ });
            }
        } else {
            LogManager::LogError("coordinate (" + std::to_string(coordi[kXIndex]) + ", "
                + std::to_string(coordi[kYIndex]) + ") is outside the computational domain in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }

        ptr_tracking_node->at(sfbitset_tmp).set_point_index.insert({0, ipoint });
        ++ipoint;
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
*/
void GeometryInfoOrigin3D::SetIndex() {
    geo_vertex_info_instance_.vec_real.resize(k0NumRealForEachVertex_);
}
int GeometryInfoOrigin3D::InitialGeometry(const DefReal dx) {
    int return_status = 0;
    this->SetIndex();
    return_status = this->GeometryInfoInterface::InitialGeometry(dx);
    return return_status;
}
/**
* @brief   function to find tracking nodes based on geometries.
* @param[in] sfbitset_aux class to manage functions for space filling code computation.
* @param[in] ptr_grid_info pointer to class store grid information.
*/
void GeometryInfoOrigin3D::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) {
    std::pair<ECriterionType, DefInt> key_tracking_grid = { ECriterionType::kGeometry, i_geo_ };
    ptr_grid_info->map_ptr_tracking_grid_info_
        .at(key_tracking_grid).get()->grid_extend_type_ = grid_extend_type_;
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
     ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset sfbitset_tmp;
    std::vector<DefReal> coordi(3, 0.);
    DefSizet ipoint = 0;
    for (const auto& iter : coordinate_origin_) {
        coordi[kXIndex] = iter.coordinate[kXIndex];
        coordi[kYIndex] = iter.coordinate[kYIndex];
        coordi[kZIndex] = iter.coordinate[kZIndex];
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(ptr_grid_info->grid_space_, coordi);
        if (ptr_grid_info->CheckIfNodeOutsideCubicDomain(3, sfbitset_tmp, sfbitset_aux)
            >= GridInfoInterface::kFlagInsideDomain_) {
            if (ptr_tracking_node->find(sfbitset_tmp) == ptr_tracking_node->end()) {
                ptr_tracking_node->insert({ sfbitset_tmp, (ptr_grid_info
                ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())->k0TrackNodeInstance_ });
            }
        } else {
            LogManager::LogError("coordinate (" + std::to_string(coordi[kXIndex]) + ", "
                + std::to_string(coordi[kYIndex]) + ", "
                + std::to_string(coordi[kZIndex]) + ") is outside the computational domain in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        ptr_tracking_node->at(sfbitset_tmp).set_point_index.insert({ 0, ipoint });
        ++ipoint;
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject