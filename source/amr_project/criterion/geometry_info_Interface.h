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
namespace grid {
class TrackingGridInfoCreatorInterface;
class GhostGridInfoCreatorInterface;
class  SFBitsetAux2D;
class  SFBitsetAux3D;
class  GridInfoInterface;
}  // end namespace grid
namespace criterion {
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
    friend class DefaultGeoShapeManager;
 public:
    // information of geometry itself
    DefReal computational_cost_ = 1.;
    DefReal decompose_factor_ = 1.;
    DefSizet i_level_ = 0;
    DefSizet i_geo_ = ~0;
    EGeometryCellType geometry_cell_type_ =
        EGeometryCellType::kUndefined;
    EGeometryEnclosedType geometry_enclosed_type_ =
        EGeometryEnclosedType::kOpen;
    std::string node_type_;

    // information stored on each vertex
    DefSizet k0NumIntForEachvertex_ = 0;
    DefSizet k0NumRealForEachvertex_ = 0;
    std::vector<GeometryvertexInfo> vec_vertexs_info_{};

    std::shared_ptr<grid::TrackingGridInfoCreatorInterface>
        ptr_tracking_grid_info_creator_ = nullptr;
    std::shared_ptr<grid::GhostGridInfoCreatorInterface>
        ptr_ghost_grid_info_creator_ = nullptr;

    // type of default geometry shape
    DefaultGeoShapeType k0DefaultGeoShapeType_ =
        DefaultGeoShapeType::kUndefined;

 protected:
    GeometryvertexInfo geo_vertex_info_instance_;
    ///< an instance for a geometry vertex with preset vector sizes
};
/**
* @class GeometryInfoInterface
* @brief class used to store information of a geometry
*/
class GeometryInfoCreatorInterface {
public:
    virtual std::shared_ptr<GeometryInfoInterface>
        CreateGeometryInfo() = 0;
};
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
* @class GeometryInfo2DInterface
* @brief class used to store information of a 2D geometry
*/
class GeometryInfo2DInterface:public GeometryInfoInterface {
public:
    // information of geometry itself
    std::array<DefReal, 2> geometry_center_{};
    std::array<DefReal, 2> flood_fill_origin_{};
    std::array<DefReal, 2> k0RealOffest_{};
    std::vector<GeometryCoordinate2D> coordinate_origin_{};

    virtual void SetIndex() = 0;
    virtual int InitialGeometry(
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo) = 0;
    virtual int UpdateGeometry(
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo) = 0;
    virtual void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) = 0;
    virtual void FindTrackingNodeNearGeo(
        const grid::SFBitsetAux2D& sfbitset_aux_2d,
        std::shared_ptr<grid::GridInfoInterface> ptr_grid_info) = 0;
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
* @class GeometryInfoInterface
* @brief class used to store information of a geometry
*/
class GeometryInfo3DInterface :public GeometryInfoInterface {
public:
    // information of geometry itself
    std::array<DefReal, 3> geometry_center_{};
    std::array<DefReal, 3> flood_fill_origin_{};
    std::array<DefReal, 3> k0RealOffest_{};
    std::vector<GeometryCoordinate3D> coordinate_origin_{};

    virtual void SetIndex() = 0;
    virtual int InitialGeometry(
        std::shared_ptr<GeometryInfo3DInterface> ptr_geo) = 0;
    virtual int UpdateGeometry(
        std::shared_ptr<GeometryInfo3DInterface> ptr_geo) = 0;
    virtual void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) = 0;
    virtual void FindTrackingNodeNearGeo(
        const grid::SFBitsetAux3D& sfbitset_aux_3d,
        std::shared_ptr<grid::GridInfoInterface> ptr_grid_info) = 0;
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

/**
* @class DefaultGeoShapeManager
* @brief class used to manage funcitons to generate default geometry shapes
* @date  2022-8-5
*/
class DefaultGeoShapeManager {
 public:
    static DefaultGeoShapeManager* GetInstance() {
        static DefaultGeoShapeManager default_geo_shape_manager_instance_;
        return &default_geo_shape_manager_instance_;
    }

    void geo_circle_initial(const DefUint dims,
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo);
    void geo_circle_update(const DefUint dims,
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo);

private:
    static DefaultGeoShapeManager* default_geo_shape_manager_instance_;
    DefaultGeoShapeManager(void) {}
    ~DefaultGeoShapeManager(void) {}
    DefaultGeoShapeManager(const DefaultGeoShapeManager&) = delete;
    DefaultGeoShapeManager& operator =
        (const DefaultGeoShapeManager&) = delete;
};
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
