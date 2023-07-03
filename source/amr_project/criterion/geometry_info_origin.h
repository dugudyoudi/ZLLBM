//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.h
* @author Zhengliang Liu
* @date  2023-4-25
* @brief  define classes to store geometry information with isolated points
*/
#ifndef ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_ORIGIN_H_
#define ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_ORIGIN_H_
#include <memory>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <utility>
#include "criterion/geometry_info_Interface.h"
#include "criterion/geometry_info_connection.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class  SFBitsetAux2D;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class  SFBitsetAux3D;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeometryInfoOrigin2D :public  Geometry2DInterface,
    public  GeometryInfoInterface {
    void SetIndex() override;
    int InitialGeometry(
        const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager) override;
    int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager) override;
    void DecomposeNHigherLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override {}
    void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface* ptr_sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex]};
    }

 public:
    GeometryInfoOrigin2D() {
        this->node_type_ = "GeometryInfoOrigin2D";
    }
    DefSizet GetNumOfGeometryPoints() const final {
            return coordinate_origin_.size();
    }
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeometryInfoOrigin3D :public  Geometry3DInterface,
    public  GeometryInfoInterface {
    void SetIndex() override;
    int InitialGeometry(
        const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager) override;
    int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager) override;
    void DecomposeNHigherLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override {}
    void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface* ptr_sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex]};
    }

 public:
    GeometryInfoOrigin3D() {
        this->node_type_ = "GeometryInfoOrigin3D";
    }
    DefSizet GetNumOfGeometryPoints() const final {
            return coordinate_origin_.size();
    }
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_ORIGIN_H_
