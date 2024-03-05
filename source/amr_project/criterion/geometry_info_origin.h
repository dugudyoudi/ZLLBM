//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
class GeometryInfoOrigin2D : public GeometryInfo2DInterface {
 public:
    void SetIndex() override;
    int InitialGeometry(const DefReal dx) override;
    void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex]};
    }
    GeometryInfoOrigin2D() {
        this->node_type_ = "GeometryInfoOrigin2D";
    }
    DefSizet GetNumOfGeometryPoints() const final {
        return coordinate_origin_.size();
    }
};
class GeometryInfoOrigin2DCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo() override {
        std::shared_ptr<GeometryInfoOrigin2D> ptr_temp =
         std::make_shared<GeometryInfoOrigin2D>();
        ptr_temp->node_type_ = "Origin2D";
        ptr_temp->geometry_cell_type_ = EGeometryCellType::kPolyLine;
        return ptr_temp;
    };
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeometryInfoOrigin3D : public GeometryInfo3DInterface {
 public:
    void SetIndex() override;
    int InitialGeometry(const DefReal dx) override;
    void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex], flood_fill_origin_[kZIndex]};
    }
    GeometryInfoOrigin3D() {
        this->node_type_ = "GeometryInfoOrigin3D";
    }
    DefSizet GetNumOfGeometryPoints() const final {
        return coordinate_origin_.size();
    }
};
class GeometryInfoOrigin3DCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo() override {
        std::shared_ptr<GeometryInfoOrigin3D> ptr_temp =
         std::make_shared<GeometryInfoOrigin3D>();
        ptr_temp->node_type_ = "Origin3D";
        ptr_temp->geometry_cell_type_ = EGeometryCellType::kPolyLine;
        return ptr_temp;
    };
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_ORIGIN_H_
