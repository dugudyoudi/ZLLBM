//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.h
* @author Zhengliang Liu
* @date  2022-9-25
* @brief  define classes to store geometry informtion with
*         connection relationship of vertices.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#define ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#include <memory>
#include <set>
#include <map>
#include <utility>
#include "criterion/geometry_info_Interface.h"
//#include "grid/grid_info_preset.h"
namespace rootproject {
namespace amrproject {
class  SFBitsetAux2D;
class  SFBitsetAux3D;
/**
* @struct GeometryConnectionSurface
* @brief structure to store surface connection information
*/
struct GeometryConnectionSurface {
    DefSizet parent_surface = ~0;
    std::vector<DefSizet> child_surface;
    //std::vector<DefSizet> edge_connection;
    std::vector<std::pair<DefSizet, DefSizet>> vertex_connection;
};
struct GeometryConnectionSurfaceLevel {
    DefSizet level_diff;
    std::vector<GeometryConnectionSurface>
        vec_surface_connection;
};
/**
* @struct GeometryConnectionEdge
* @brief structure to store edge connection information
*/
struct GeometryConnectionEdge {
    //DefUint status_added = 2;
    //std::array<std::pair<DefSizet, DefSizet>, 2> vertex_connection;
    std::set<DefSizet> set_index_surfaces;
};

struct GeometryConnectionEdgeLevel {
    DefSizet level_diff;
    // the first pair<DefSizet, DefSizet> is the vertex whose index is larger
    std::map<std::pair<std::pair<DefSizet, DefSizet>,
        std::pair<DefSizet, DefSizet>>, GeometryConnectionEdge>
        map_edge_connection;
};
/**
* @struct GeometryConnectionCoordinate
* @brief structure to store vertex formation
*/
struct GeometryConnectionCoordinate {
    //bool bool_decomposed = false;
    std::map<DefSizet,
        std::set<std::pair<DefSizet, DefSizet>>> map_linked_vertices_level;
    std::array<std::pair<DefSizet, DefSizet>, 2> parent_vertices;
    std::set<std::pair<DefSizet, DefSizet>> child_vertices;
    GeometryVertexInfo vertex_info;
    std::vector<DefReal> coordinates;
    std::map<DefSizet, DefSFBitset> map_bitset_ref;
    DefSizet highest_grid_level = 0;
};
/**
* @struct GeometryConnectionCoordinateLevel
* @brief structure to store all vertex formation for a level
*/
struct GeometryConnectionCoordinateLevel {
    std::vector<GeometryConnectionCoordinate> vec_vertex_coordinate;
};
/**
* @class GeometryConnectionInterface
* @brief class for geometries described by connection relations
*/
class GeometryConnectionInterface : virtual public GeometryInfoInterface {
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
    // vertex_given_level_ remain the same all the time
    bool bool_change_vertices_from_geo_ = false;
    bool bool_vertex_info_stored_for_connection_ = false;
    DefSizet k0NumIntForEachVertex_ = 0;
    DefSizet k0NumRealForEachVertex_ = 0;
    DefSizet k0NumEdgeForSurface_ = 3;

    std::vector<std::vector<DefSizet>> connection_relation_;

    std::vector<GeometryConnectionCoordinateLevel>
        vertex_given_level_{};
    std::vector<GeometryConnectionEdgeLevel>
        connection_edge_given_level_{};
    std::vector<GeometryConnectionSurfaceLevel>
        connection_surface_given_level_{};
    std::vector<std::set<std::pair<DefSizet, DefSizet>>>
        connection_vertex_given_level_{};
        
    virtual void InitialCoordinateGivenLevel(
        std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max) = 0;
    virtual DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1) = 0;
    virtual void ComputeMidCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1,
        std::vector<DefReal>* const ptr_coordinates) = 0;
    virtual ~GeometryConnectionInterface() {}

    void SetupConnectionParameters(EGeometryCellType cell_type);
    void InitialConnection(std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max);
    void MergeEdgeOnce(const DefSizet i_input_level, const DefReal ds_min,
        const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
        const std::set<std::pair<std::pair<DefSizet, DefSizet>,
        std::pair<DefSizet, DefSizet>>>& edge_for_merge,
        std::set<std::pair<std::pair<DefSizet, DefSizet>,
        std::pair<DefSizet, DefSizet>>>* const  ptr_surface_remain_for_merge,
        DefMap<DefUint>* const ptr_sfbitset_ref_removed);
    void BisectEdgeOnce(const DefSizet i_input_level, const DefReal ds_max,
        const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
        const std::set<std::pair<std::pair<DefSizet, DefSizet>,
        std::pair<DefSizet, DefSizet>>>& edge_for_bisect,
        std::set<std::pair<std::pair<DefSizet, DefSizet>,
        std::pair<DefSizet, DefSizet>>>* const  ptr_surface_remain_for_biset,
        DefMap<DefUint>* const ptr_sfbitset_ref_added);
    
    // virtual function from GeometryInfoInterface
    virtual void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface* ptr_sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) override final;

 protected:
    GeometryConnectionCoordinate vertex_instance_;
    void RemoveVertex(const DefSizet i_input_level,
        const std::vector<DefReal>& grid_space,
        const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
        const std::set<std::pair<DefSizet, DefSizet>>& set_vertex_remove,
        DefMap<DefUint>* const ptr_sfbitset_ref_removed);
    void AddNewLinkage(const DefSizet i_input_level,
        const std::pair<DefSizet, DefSizet>& vertex_new,
        const std::pair<DefSizet, DefSizet>& vertex_origin);
    void FindSurfaceForReconstruction(const DefSizet i_input_level,
        const std::set<DefSizet>& surface_process,
        const std::set<std::pair<DefSizet, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_reconstruct);
    void ReconstructSurfaceBasedOnExistingVertex(const DefSizet i_input_level,
        const std::set<DefSizet>& surface_reconstruct,
        const std::set<std::pair<DefSizet, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_remove);
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeometryInfoConnection2D:public Geometry2DInterface,
    public GeometryConnectionInterface {
 public:
    DefSizet k0UxIndex_ = 0, k0UyIndex_ = 0;
    DefSizet k0FxIndex_ = 0, k0FyIndex_ = 0;

    // virtual functions for Geometry2DInterface
    void SetIndex() override;
    int InitialGeometry(const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_managerr) override;
    int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager) override;
    void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override;

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel(
        std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max) override;
    DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1) override;
    void ComputeMidCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1,
        std::vector<DefReal>* const ptr_coordinates) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const override final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex] };
    }
    DefSizet GetNumOfGeometryPoints() const override final {
        return coordinate_origin_.size();
    }
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
class GeometryInfoConnection3D :public Geometry3DInterface,
    public GeometryConnectionInterface {
 public:
    DefSizet k0UxIndex_ = 0, k0UyIndex_ = 0, k0UzIndex_ = 0;
    DefSizet k0FxIndex_ = 0, k0FyIndex_ = 0, k0FzIndex_ = 0;

    // virtual functions for Geometry3DInterface
    void SetIndex() override;
    int InitialGeometry(const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager) override;
    int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager) override;
    void DecomposeNHigerLevel(const DefSizet i_level_grid,
        const DefReal decompose_length,
        const std::unordered_map<DefSizet, bool>& map_indices_base,
        std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain)
        override;

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel(
        std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max) override;
    DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1) override;
    void ComputeMidCoordinates(
        const std::pair<DefSizet, DefSizet>& vertex0,
        const std::pair<DefSizet, DefSizet>& vertex1,
        std::vector<DefReal>* const ptr_coordinates) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const override final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex],
        flood_fill_origin_[kZIndex] };
    }
    DefSizet GetNumOfGeometryPoints() const override final {
        return coordinate_origin_.size();
    }
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
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
