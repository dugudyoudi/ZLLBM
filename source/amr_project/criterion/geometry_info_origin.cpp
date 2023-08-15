//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.cpp
* @author Zhengliang Liu
* @brief functions for geometries with connection relations, e.g. stl files
* @date  2023-4-25
*/
#include <vector>
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
int GeometryInfoOrigin2D::InitialGeometry(const DefReal dx,
    const DefaultGeoShapeType shape_type,
    const DefaultGeoManager& default_geo_manager) {
    int return_status = 0;
    this->SetIndex();
    this->k0DefaultGeoShapeType_ = shape_type;
    return_status = this->GeometryInfo2DInterface::InitialGeometry(
        dx, shape_type, default_geo_manager);
    return return_status;
}
int GeometryInfoOrigin2D::UpdateGeometry(
    const DefaultGeoManager& default_geo_manager) {
    return 0;
}
void GeometryInfoOrigin2D::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface* ptr_sfbitset_aux,
    GridInfoInterface* const ptr_grid_info) {
    std::pair<ECriterionType, DefAmrIndexUint> key_tracking_grid = { ECriterionType::kGeometry, i_geo_ };
    if (ptr_grid_info->map_ptr_tracking_grid_info_.find(key_tracking_grid)
     == ptr_grid_info->map_ptr_tracking_grid_info_.end()) {
        ptr_grid_info->map_ptr_tracking_grid_info_.insert(
         { key_tracking_grid, ptr_tracking_grid_info_creator_->CreateTrackingGridInfo() });
    }
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
     ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset bitset_temp;
    std::vector<DefReal> coordi(2, 0.);
    DefSizet ipoint = 0;
    for (const auto& iter : coordinate_origin_) {
        coordi[kXIndex] = iter.coordinate[kXIndex];
        coordi[kYIndex] = iter.coordinate[kYIndex];
        bitset_temp = ptr_sfbitset_aux->SFBitsetEncodingCoordi(ptr_grid_info->grid_space_, coordi);
        if (ptr_tracking_node->find(bitset_temp) == ptr_tracking_node->end()) {
            ptr_tracking_node->insert({ bitset_temp, (ptr_grid_info
             ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())->k0TrackNodeInstance_ });
        }
        ptr_tracking_node->at(bitset_temp).set_point_index.insert({0, ipoint });
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
int GeometryInfoOrigin3D::InitialGeometry(const DefReal dx,
    const DefaultGeoShapeType shape_type,
    const DefaultGeoManager& default_geo_manager) {
    int return_status = 0;
    this->SetIndex();
    this->k0DefaultGeoShapeType_ = shape_type;
    return_status = this->GeometryInfo3DInterface::InitialGeometry(
        dx, shape_type, default_geo_manager);
    return return_status;
}
int GeometryInfoOrigin3D::UpdateGeometry(
    const DefaultGeoManager& default_geo_manager) {
    return 0;
}
void GeometryInfoOrigin3D::FindTrackingNodeBasedOnGeo(
    const SFBitsetAuxInterface* ptr_sfbitset_aux,
    GridInfoInterface* const ptr_grid_info) {
    std::pair<ECriterionType, DefAmrIndexUint> key_tracking_grid =
    { ECriterionType::kGeometry, i_geo_ };
    if (ptr_grid_info->map_ptr_tracking_grid_info_.find(key_tracking_grid)
     == ptr_grid_info->map_ptr_tracking_grid_info_.end()) {
        ptr_grid_info->map_ptr_tracking_grid_info_
         .insert({ key_tracking_grid,
         ptr_tracking_grid_info_creator_->CreateTrackingGridInfo() });
    }
    DefMap<TrackingNode>* ptr_tracking_node = &(ptr_grid_info
     ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get()->map_tracking_node_);

    DefSFBitset bitset_temp;
    std::vector<DefReal> coordi(3, 0.);
    DefSizet ipoint = 0;
    for (const auto& iter : coordinate_origin_) {
        coordi[kXIndex] = iter.coordinate[kXIndex];
        coordi[kYIndex] = iter.coordinate[kYIndex];
        coordi[kZIndex] = iter.coordinate[kZIndex];
        bitset_temp = ptr_sfbitset_aux->SFBitsetEncodingCoordi(
            ptr_grid_info->grid_space_, coordi);
        if (ptr_tracking_node->find(bitset_temp) == ptr_tracking_node->end()) {
            ptr_tracking_node->insert({ bitset_temp, (ptr_grid_info
             ->map_ptr_tracking_grid_info_.at(key_tracking_grid).get())
             ->k0TrackNodeInstance_ });
        }
        ptr_tracking_node->at(bitset_temp).set_point_index
            .insert({ 0, ipoint });
        ++ipoint;
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject