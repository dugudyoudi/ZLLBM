//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define classes to store geometry information.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#define ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#include <vector>
#include <array>
#include <memory>
#include <string>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
namespace rootproject {
namespace amrproject {
class TrackingGridInfoCreatorInterface;
class GhostGridInfoCreatorInterface;
class SFBitsetAuxInterface;
class GridInfoInterface;
/**
* @struct GeometryvertexInfo
* @brief struct used to store information of a geometry vertex
*/
struct GeometryvertexInfo {
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
    DefReal computational_cost_ = 1.;
    DefReal decompose_factor_ = 1.;
    DefSizet i_level_ = 0;
    DefSizet i_geo_ = ~0;
    EGeometryCellType geometry_cell_type_ =
        EGeometryCellType::kUndefined;
    EGridExtendType grid_extend_type_ =
        EGridExtendType::kSameInAllDirections;
    std::string node_type_;

    // information stored on each vertex
    DefSizet k0NumIntForEachvertex_ = 0;
    DefSizet k0NumRealForEachvertex_ = 0;
    std::vector<GeometryvertexInfo> vec_vertices_info_{};

    std::shared_ptr<TrackingGridInfoCreatorInterface>
        ptr_tracking_grid_info_creator_ = nullptr;
    std::shared_ptr<GhostGridInfoCreatorInterface>
        ptr_ghost_grid_info_creator_ = nullptr;

    // type of default geometry shape
    DefaultGeoShapeType k0DefaultGeoShapeType_ =
        DefaultGeoShapeType::kUndefined;

    // number of extended layer based on geometry
    std::vector<DefLUint>
        k0XIntExtendPositive_, k0XIntExtendNegative_,
        k0YIntExtendPositive_, k0YIntExtendNegative_,
        k0ZIntExtendPositive_, k0ZIntExtendNegative_;
        ///< number of extened layers
   
    /* number of layer extended inside the geometry
     at (i_level - 1) refienemnt level*/
    std::vector<DefLUint>  k0IntInnerExtend_;
    ///< number of extened layers inside the geometry

    virtual void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface* ptr_sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) = 0;
    virtual std::vector<DefReal> GetFloodFillOriginArrAsVec() const = 0;
    virtual DefSizet GetNumOfGeometryPoints() const = 0;
    virtual ~GeometryInfoInterface() {}
 protected:
    GeometryvertexInfo geo_vertex_info_instance_;
    ///< instance for a geometry vertex with preset vector sizes
};
/**
* @class GeometryInfoCreatorInterface
* @brief interface class used to create geometry instance
*/
class GeometryInfoCreatorInterface {
public:
    virtual std::shared_ptr<GeometryInfoInterface>
        CreateGeometryInfo() = 0;
};
class DefaultGeoManager;
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @struct GeometryCoordinate2D
* @brief struct used to store information of a geometry vertex
*/
struct GeometryCoordinate2D {
    std::array<DefReal, 2> coordinate{};

    //GeometryCoordinate2D& operator=(const GeometryCoordinate2D& coordi_r) {
    //    this->coordinate.at(0) = coordi_r.coordinate.at(0);
    //    this->coordinate.at(1) = coordi_r.coordinate.at(1);
    //    return *this;
    //}
    //GeometryCoordinate2D operator+(const GeometryCoordinate2D& coordi_r) {
    //    GeometryCoordinate2D coordi_return;
    //    coordi_return.coordinate.at(0) = this->coordinate.at(0)
    //        + coordi_r.coordinate.at(0);
    //    coordi_return.coordinate.at(1) = this->coordinate.at(1)
    //        + coordi_r.coordinate.at(1);
    //    return coordi_return;
    //}
    //GeometryCoordinate2D operator/(DefReal real_r) {
    //    GeometryCoordinate2D coordi_return;
    //    coordi_return.coordinate.at(0) = this->coordinate.at(0) / real_r;
    //    coordi_return.coordinate.at(1) = this->coordinate.at(1) / real_r;
    //}
};
/**
* @class Geometry2DInterface
* @brief interface class used to store 2D information of a geometry
*/
class Geometry2DInterface{
public:
    // information of geometry itself
    std::array<DefReal, 2> geometry_center_{};
    std::array<DefReal, 2> flood_fill_origin_{};
    std::array<DefReal, 2> k0RealOffest_{};
    std::vector<GeometryCoordinate2D> coordinate_origin_{};

    virtual int InitialGeometry(const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager);
    virtual int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager);

    virtual void SetIndex() = 0;
    virtual void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) = 0;
    virtual ~Geometry2DInterface() {}
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @struct GeometryCoordinate3D
* @brief struct used to store information of a geometry vertex
*/
struct GeometryCoordinate3D {
    std::array<DefReal, 3> coordinate{};
};
/**
* @class Geometry3DInterface
* @brief interface class used to store 3D information of a geometry
*/
class Geometry3DInterface {
public:
    // information of geometry itself
    std::array<DefReal, 3> geometry_center_{};
    std::array<DefReal, 3> flood_fill_origin_{};
    std::array<DefReal, 3> k0RealOffest_{};
    std::vector<GeometryCoordinate3D> coordinate_origin_{};

    virtual void SetIndex() = 0;
    virtual int InitialGeometry(const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager);
    virtual int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager);
    virtual void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) = 0;
    virtual ~Geometry3DInterface() {}
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
class DefaultGeoManager {
public:
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    void circle_initial(Geometry2DInterface* const ptr_geo) const;
    void circle_update(DefReal sum_t,
        Geometry2DInterface* const ptr_geo) const;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
    void cube_initial(const DefReal dx,
        Geometry3DInterface* const ptr_geo) const;
    void cube_update(const DefReal sum_t,
        Geometry3DInterface* const ptr_geo) const;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
