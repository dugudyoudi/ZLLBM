//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.h
* @author Zhengliang Liu
* @date  2022-9-25
* @brief  define classes to store geometry informtion with
*         connection relationship of points.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#define ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#include <memory>
#include <set>
#include <utility>
#include "criterion/geometry_info_Interface.h"
//#include "grid/grid_info_preset.h"
namespace rootproject {
namespace amrproject {
namespace criterion {
/**
* @struct GeometryConnectionCoordinateInterface
* @brief structure to store index of vertex
*/
struct GeometryCoordinateIndex {
    DefSizet layer_level = 0; /**< difference between the level of edge
                           and the level of coordinate*/
    DefSizet vertex_index;
    GeometryCoordinateIndex& operator=(
        const GeometryCoordinateIndex& coordi_r) {
        this->layer_level = coordi_r.layer_level;
        this->vertex_index = coordi_r.vertex_index;
        return *this;
    }
};
/**
* @struct GeometryConnectionSurface
* @brief structure to store surface connection information
*/
struct GeometryConnectionSurface {
    DefSizet parent_surface;
    std::vector<DefSizet> child_surface;
    std::vector<DefSizet> edge_connection;
    std::vector<GeometryCoordinateIndex> vertex_connectin;
};
struct GeometryConnectionSurfaceLevel {
    DefSizet i_level;
    std::vector<GeometryConnectionSurface>
        vec_surface_connection;
};
/**
* @struct GeometryConnectionEdge
* @brief structure to store edge connection information
*/
struct GeometryConnectionEdge {
    //DefUint status_added = 2;
    std::array<GeometryCoordinateIndex, 2> vertex_connection;
    std::vector<DefSizet> vec_index_surfaces;
};
struct GeometryConnectionEdgeLevel {
    DefSizet i_level;
    std::vector<GeometryConnectionEdge>
        vec_edge_connection;
};
/**
* @struct GeometryConnectionCoordinateInterface
* @brief structure to store vertex formation
*/
struct GeometryConnectionCoordinateInterface {
    //bool bool_decomposed = false;
    std::vector<std::pair<DefSizet, DefSizet>> vec_index_edges;
    std::array<GeometryCoordinateIndex, 2> parent_vertices;
    std::vector<GeometryCoordinateIndex> vec_child_vertices;
    GeometryPointInfo point_info;
};
struct GeometryConnectionCoordinate2D :public GeometryCoordinate2D,
    public GeometryConnectionCoordinateInterface {
};
struct GeometryConnectionCoordinate3D :public GeometryCoordinate3D,
    public GeometryConnectionCoordinateInterface {
};
struct GeometryConnectionCoordinate2DLevel {
    std::vector<GeometryConnectionCoordinate2D>
        vec_vertex_cooridinate;
};
struct GeometryConnectionCoordinate3DLevel {
    std::vector<GeometryConnectionCoordinate3D>
        vec_vertex_cooridinate;
};
class GeometryConnectionInterface  {
    // index of parameters in GeometryCoordinate::vec_real
    // k0 indicates those parameters are optional
 public:
    bool bool_vec_velocities_ = true;
    bool bool_vec_forces_ = true;

    // if true, the last and first elements in the connection relation 
    // (connection_relation_) consist of an edge
    bool bool_periodic_connection_ = false;

    // if false, vertices from coordinate_origin_ will not change,
    // i.e. the first coordinate_origin_.size() elements in
    // coordinate_given_level_ remain the same all the time
    bool bool_change_vertices_from_geo_ = false;
    bool bool_vertex_info_ = false;
    DefSizet k0NumIntForEachVertex_ = 0;
    DefSizet k0NumRealForEachVertex_ = 0;
    DefSizet k0NumEdgeForSurface_ = 3;

    std::vector<std::vector<DefSizet>> connection_relation_;

    std::vector<GeometryConnectionEdgeLevel>
        connection_edge_given_level_{};
    std::vector<GeometryConnectionSurfaceLevel>
        connection_surface_given_level_{};
    std::vector<std::set<std::pair<DefSizet, DefSizet>>>
        connection_index_given_level_{};
        
    virtual void InitialCoordinateGivenLevel() {
        printf_s("Initialize 2D or 3D coordinates"
            " with dimension specified function.");
    }
    virtual void InitialConnectionAddEdge(const DefSizet index_small,
        const DefSizet index_large, const DefSizet i_surface) {
        printf_s("Add edges and connection relation during initialization"
            " with dimension specified function.");
    }

    void SetupConnectionParameters(EGeometryCellType cell_type);
    void InitialConnectionGivenLevel(DefSizet i_level);
    void DecomposeOneHigerLevel(
        const DefSizet i_level, const DefReal decompose_length);
 protected:
    GeometryPointInfo connection_vertex_info_instance_;
    GeometryConnectionSurface connection_surface_instance_;
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeometryInfoConnection2D:public GeometryInfo2DInterface,
    public GeometryConnectionInterface {
 public:
    DefSizet k0UxIndex_ = 0, k0UyIndex_ = 0;
    DefSizet k0FxIndex_ = 0, k0FyIndex_ = 0;
    std::vector<GeometryConnectionCoordinate2DLevel>
        coordinate_given_level_{};

    // virtual functions for GeometryInfoInterface
    void SetIndex() override;
    int InitialGeometry(
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo) override;
    int UpdateGeometry(
        std::shared_ptr<GeometryInfo2DInterface> ptr_geo) override;
    void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override;
    void DecomposeOneHigerLevel(
        const DefSizet i_level, const DefReal decompose_length);

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel() override;
    void InitialConnectionAddEdge(const DefSizet index_small,
        const DefSizet index_large, const DefSizet i_surface) override;
    void BisectinOnceLine(const DefSizet i_input_level, const DefReal ds_max,
        const std::vector<DefSizet>& vec_edge_for_bisect,
        std::vector<DefSizet>* ptr_surface_remain_for_bisect);
    void MergeOnceLine(const DefSizet i_input_level, const DefReal ds_min,
        const std::vector<DefSizet>& vec_edge_for_bisect,
        std::vector<DefSizet>* ptr_surface_remain_for_bisect);
    void MergeOnceTriangle(const DefSizet i_input_level, const DefReal ds_min,
        const std::set<DefSizet>& vec_edge_for_bisect,
        std::set<DefSizet>* ptr_surface_remain_for_merge);
private:
    void DeleteChildRelation(const DefSizet layer_level, const DefSizet index);
};
class GeometryInfoConnection2DCreator :public GeometryInfoCreatorInterface {
public:
    std::shared_ptr<GeometryInfoInterface>
        CreateGeometryInfo() override {
        std::shared_ptr<GeometryInfoConnection2D> ptr_temp =
            std::make_shared<GeometryInfoConnection2D>();
        ptr_temp->node_type_ = "Connection2D";
        ptr_temp->geometry_cell_type_ = EGeometryCellType::kPolyLine;
        return ptr_temp;
    };
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeometryInfoConnection3D :public GeometryInfo3DInterface,
    GeometryConnectionInterface {
 public:
    DefSizet k0UxIndex_ = 0, k0UyIndex_ = 0, k0UzIndex_ = 0;
    DefSizet k0FxIndex_ = 0, k0FyIndex_ = 0, k0FzIndex_ = 0;
    std::vector<GeometryConnectionCoordinate3DLevel>
        coordinate_given_level_{};

    void SetIndex() override;
    int InitialGeometry(
        std::shared_ptr<GeometryInfo3DInterface> ptr_geo) override;
    int UpdateGeometry(
        std::shared_ptr<GeometryInfo3DInterface> ptr_geo) override;
    void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override;

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel() override;
    void InitialConnectionAddEdge(const DefSizet index_small,
        const DefSizet index_large, const DefSizet i_surface) override {
    };
};
class GeometryInfoConnection3DCreator :public GeometryInfoCreatorInterface {
public:
    std::shared_ptr<GeometryInfoInterface>
        CreateGeometryInfo() override {
        std::shared_ptr<GeometryInfoConnection3D> ptr_temp =
            std::make_shared<GeometryInfoConnection3D>();
        ptr_temp->node_type_ = "Connection3D";
        ptr_temp->geometry_cell_type_ = EGeometryCellType::kPolyLine;
        return ptr_temp;
    };
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
