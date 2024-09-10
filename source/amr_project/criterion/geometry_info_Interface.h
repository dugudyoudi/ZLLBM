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
#include "criterion/geometry_coordi.h"
#include "criterion/geometry_default_shape.h"
namespace rootproject {
namespace amrproject {
class TrackingGridInfoCreatorInterface;
class GhostGridInfoCreatorInterface;
class SFBitsetAuxInterface;
class GridInfoInterface;
class GeoShapeInterface;
/**
* @struct GeometryVertexInfo
* @brief struct used to store information of a geometry vertex
*/
struct GeometryVertexInfo {
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
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
    DefInt k0NumIntForEachVertex_ = 0;
    DefInt k0NumRealForEachVertex_ = 0;
    std::vector<GeometryVertexInfo> vec_vertices_info_{};

    TrackingGridInfoCreatorInterface* ptr_tracking_grid_info_creator_ = nullptr;
    GhostGridInfoCreatorInterface* ptr_ghost_grid_info_creator_ = nullptr;


    // number of extended layer based on geometry
    std::vector<DefInt> k0XIntExtendPositive_, k0XIntExtendNegative_,
     k0YIntExtendPositive_, k0YIntExtendNegative_, k0ZIntExtendPositive_, k0ZIntExtendNegative_;
     ///< number of extened layers

    /* number of layer extended inside the geometry
     at (i_level - 1) refinement level*/
    std::vector<DefInt>  k0IntInnerExtend_;
    ///< number of extened layers inside the geometry

    virtual void SetIndex() = 0;
    virtual bool SetOffset(const std::vector<DefReal>& vec_offset) = 0;
    std::unique_ptr<GeoShapeInterface> ptr_geo_shape_;
    virtual int InitialGeometry(const DefReal dx);
    virtual int UpdateGeometry(const DefReal sum_t);
    virtual void FindTrackingNodeBasedOnGeo(
     const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) = 0;
    virtual std::vector<DefReal> GetFloodFillOriginArrAsVec() const = 0;
    virtual DefSizet GetNumOfGeometryPoints() const = 0;
    virtual ~GeometryInfoInterface() {}

 protected:
    GeometryVertexInfo geo_vertex_info_instance_;
    ///< instance for a geometry vertex with preset vector sizes
};
/**
* @class GeometryInfoCreatorInterface
* @brief interface class used to create geometry instance
*/
class GeometryInfoCreatorInterface {
 public:
    virtual std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo() = 0;
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @class GeometryInfo2DInterface
* @brief interface class used to store 2D information of a geometry
*/
class GeometryInfo2DInterface: virtual public GeometryInfoInterface {
 public:
    // information of geometry itself
    std::array<DefReal, 2> geometry_center_{};
    std::array<DefReal, 2> flood_fill_origin_{};
    std::array<DefReal, 2> k0RealMin_{};
    std::vector<GeometryCoordinate2D> coordinate_origin_{};

    bool SetOffset(const std::vector<DefReal>& vec_offset) final {
        if (vec_offset.size() != 2) {
            return false;
        }
        k0RealMin_[kXIndex] = vec_offset.at(kXIndex);
        k0RealMin_[kYIndex] = vec_offset.at(kYIndex);
        return true;
    }

    virtual ~GeometryInfo2DInterface() {}
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @class GeometryInfo3DInterface
* @brief interface class used to store 3D information of a geometry
*/
class GeometryInfo3DInterface: virtual public GeometryInfoInterface {
 public:
    // information of geometry itself
    std::array<DefReal, 3> geometry_center_{};
    std::array<DefReal, 3> flood_fill_origin_{};
    std::array<DefReal, 3> k0RealMin_{};
    std::vector<GeometryCoordinate3D> coordinate_origin_{};

    bool SetOffset(const std::vector<DefReal>& vec_offset) final {
        if (vec_offset.size() != 3) {
            return false;
        }
        k0RealMin_[kXIndex] = vec_offset.at(kXIndex);
        k0RealMin_[kYIndex] = vec_offset.at(kYIndex);
        k0RealMin_[kZIndex] = vec_offset.at(kZIndex);
        return true;
    }

    virtual ~GeometryInfo3DInterface() {}
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
