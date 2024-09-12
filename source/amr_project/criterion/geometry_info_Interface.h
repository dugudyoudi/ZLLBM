//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define classes to store geometry information.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#define SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#include <vector>
#include <array>
#include <memory>
#include <unordered_map>
#include <string>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
#include "criterion/geometry_default_shape.h"
namespace rootproject {
namespace amrproject {
class TrackingGridInfoCreatorInterface;
class GhostGridInfoCreatorInterface;
class SFBitsetAuxInterface;
class GridInfoInterface;
class GeoShapeInterface;
/**
* @struct GeometryVertex
* @brief struct used to store information of a geometry vertex
*/
struct GeometryVertex {
    std::array<DefReal, 3> coordinate{};
};
/**
* @class GeometryInfoInterface
* @brief class used to store information of a geometry
*/
class GeometryInfoInterface {
 public:
    // information of geometry itself
    DefInt computational_cost_ = 1;
    DefReal decompose_factor_ = 1.;
    DefInt i_level_ = 0;
    DefInt i_geo_ = ~0;
    EGeometryCellType geometry_cell_type_ = EGeometryCellType::kUndefined;
    EGridExtendType grid_extend_type_ = EGridExtendType::kSameInAllDirections;
    EGeometryStatus geometry_status_ = EGeometryStatus::kVirtual;
    std::string node_type_;

    // information stored on each vertex
    std::vector<std::unique_ptr<GeometryVertex>> vec_vertices_{};
    virtual DefSizet GetNumOfGeometryPoints() const {
        return vec_vertices_.size();
    }

    TrackingGridInfoCreatorInterface* ptr_tracking_grid_info_creator_ = nullptr;
    GhostGridInfoCreatorInterface* ptr_ghost_grid_info_creator_ = nullptr;

    // number of extended layer based on geometry
    std::vector<DefInt> k0XIntExtendPositive_, k0XIntExtendNegative_, k0YIntExtendPositive_,
        k0YIntExtendNegative_, k0ZIntExtendPositive_, k0ZIntExtendNegative_;
        ///< number of extened layers

    /* number of layer extended inside the geometry
     at (i_level - 1) refinement level*/
    std::vector<DefInt> k0IntInnerExtend_;
    ///< number of extened layers inside the geometry

    void SetOffset(const std::array<DefReal, 3>& array_offset) {
        k0RealMin_ = array_offset;
    }
    std::unique_ptr<GeoShapeInterface> ptr_geo_shape_;
    virtual int InitialGeometry(const DefReal dx);
    virtual int UpdateGeometry(const DefReal sum_t);
    virtual void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) = 0;

    std::vector<DefReal> GetFloodFillOriginArrAsVec() const {
        if (k0GeoDim_ == 2) {
            return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex]};
        } else {
            return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex], flood_fill_origin_[kZIndex]};
        }
    }

    std::array<DefReal, 3> geometry_center_{};
    std::array<DefReal, 3> flood_fill_origin_{};
    std::array<DefReal, 3> k0RealMin_{};
    virtual std::unique_ptr<GeometryVertex> GeoVertexCreator() {
        return std::make_unique<GeometryVertex>();
    }

    explicit GeometryInfoInterface(const DefInt dims) : k0GeoDim_(dims) {
        this->node_type_ = "GeometryInfoInterface";
        this->geometry_cell_type_ = EGeometryCellType::kPolyLine;
    }
    virtual ~GeometryInfoInterface() {}

 protected:
    DefInt k0GeoDim_ = 0;

 private:
    GeometryInfoInterface();
};
/**
* @class GeometryInfoCreatorInterface
* @brief interface class used to create geometry instance
*/
class GeometryInfoCreatorInterface {
 public:
    virtual std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo(const DefInt dims) = 0;
};

}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
