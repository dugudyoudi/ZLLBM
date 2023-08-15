//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.h
* @author Zhengliang Liu
* @date  2022-9-25
* @brief  define classes to store geometry information with
*         connection relationship of vertices.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#define ROOTPROJECT_AMR_PROJECT_GEOMETRY_GEOMETRY_INFO_CONNECTION_H_
#include <memory>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <utility>
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class  SFBitsetAux2D;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class  SFBitsetAux3D;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
/**
* @struct GeometryConnectionSurface
* @brief structure to store surface connection information
*/
struct GeometryConnectionSurface {
    DefSizet parent_surface = ~0;
    std::vector<DefSizet> child_surface;
    std::vector<std::pair<DefAmrIndexUint, DefSizet>> vertex_connection;
};
struct GeometryConnectionSurfaceLevel {
    DefAmrIndexUint level_diff;
    std::vector<GeometryConnectionSurface>
        vec_surface_connection;
};
/**
* @struct GeometryConnectionEdge
* @brief structure to store edge connection information
*/
struct GeometryConnectionEdge {
    std::set<DefSizet> set_index_surfaces;
};

struct GeometryConnectionEdgeLevel {
    DefAmrIndexUint level_diff;
    // the first pair<DefAmrIndexUint, DefSizet> is the vertex whose index is larger
    std::map<std::pair<std::pair<DefAmrIndexUint, DefSizet>,
     std::pair<DefAmrIndexUint, DefSizet>>, GeometryConnectionEdge> map_edge_connection;
};
/**
* @struct GeometryConnectionCoordinate
* @brief structure to store vertex formation
*/
struct GeometryConnectionCoordinate {
    std::map<DefAmrIndexUint, std::set<std::pair<DefAmrIndexUint, DefSizet>>> map_linked_vertices_level;
    ///< indices of vertices at current or different levels linked to this vertex
    std::array<std::pair<DefAmrIndexUint, DefSizet>, 2> parent_vertices;
    std::set<std::pair<DefAmrIndexUint, DefSizet>> child_vertices;
    GeometryVertexInfo vertex_info;
    std::vector<DefReal> coordinates;
    std::map<DefAmrIndexUint, DefSFBitset> map_bitset_ref;  ///< spacing filing code of vertices
    DefAmrIndexUint highest_grid_level = 0;
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
class GeometryConnectionInterface {
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
    DefAmrIndexUint k0NumEdgeForSurface_ = 3;

    std::vector<std::vector<DefSizet>> connection_relation_;

    std::vector<GeometryConnectionCoordinateLevel>
        vertex_given_level_{};   ///< vertices at the given level (the ith element)
    std::vector<GeometryConnectionEdgeLevel>
        connection_edge_given_level_{};  ///< edges at the given level (the ith element)
    std::vector<GeometryConnectionSurfaceLevel>
        connection_surface_given_level_{};  ///< surfaces at the given level (ith element)
    // connection_vertex_given_level_ records the vertices at a given level since vertices may exist simultaneously
    // at lower and higher levels, while vertex_given_level_ only store vertices at the lowest levels to
    // to reduce memory cost
    std::vector<std::set<std::pair<DefAmrIndexUint, DefSizet>>>
        connection_vertex_given_level_{};
    ///<  vertices at the current and higher geometry levels exist simultaneously at the given level (the ith element)

    virtual void InitialCoordinateGivenLevel(
     std::vector<DefReal>* const ptr_coordi_min, std::vector<DefReal>* const ptr_coordi_max) = 0;
    virtual DefReal ComputeDistanceFromCoordinates(
     const std::pair<DefAmrIndexUint, DefSizet>& vertex0, const std::pair<DefAmrIndexUint, DefSizet>& vertex1) = 0;
    virtual void ComputeMidCoordinates(const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
     const std::pair<DefAmrIndexUint, DefSizet>& vertex1, std::vector<DefReal>* const ptr_coordinates) = 0;
    virtual ~GeometryConnectionInterface() {}

    void SetupConnectionParameters(EGeometryCellType cell_type);
    void InitialConnection(std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max);
    void MergeEdgeOnce(const DefAmrIndexUint i_level, const DefAmrIndexUint i_input_level, const DefReal ds_min,
     const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
     const std::set<std::pair<std::pair<DefAmrIndexUint, DefSizet>,
      std::pair<DefAmrIndexUint, DefSizet>>>& edge_for_merge,
     std::set<std::pair<std::pair<DefAmrIndexUint, DefSizet>, std::pair<DefAmrIndexUint, DefSizet>>>*
     const ptr_edhe_remain_for_merge, DefMap<DefAmrUint>* const ptr_sfbitset_ref_removed);
    void BisectEdgeOnce(const DefAmrIndexUint i_level, const DefAmrIndexUint i_input_level, const DefReal ds_max,
     const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
     const std::set<std::pair<std::pair<DefAmrIndexUint, DefSizet>,
      std::pair<DefAmrIndexUint, DefSizet>>>& edge_for_bisect,
     std::set<std::pair<std::pair<DefAmrIndexUint, DefSizet>, std::pair<DefAmrIndexUint, DefSizet>>>*
     const ptr_surface_remain_for_bisect, DefMap<DefAmrUint>* const ptr_sfbitset_ref_added);

    void FindTrackingNodeBasedOnGeo(DefAmrIndexUint i_geo, DefAmrIndexUint i_level,
     const EGridExtendType grid_extend_type, const TrackingGridInfoCreatorInterface& ptr_tracking_creator,
     const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info);

 protected:
    GeometryConnectionCoordinate vertex_instance_;
    void RemoveVertex(const DefAmrIndexUint i_input_level,
        const std::vector<DefReal>& grid_space,
        const std::shared_ptr<SFBitsetAuxInterface> ptr_sfbitset_aux,
        const std::set<std::pair<DefAmrIndexUint, DefSizet>>& set_vertex_remove,
        DefMap<DefAmrUint>* const ptr_sfbitset_ref_removed);
    void AddNewLinkage(const DefAmrIndexUint i_input_level,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex_new,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex_origin);
    void FindSurfaceForReconstruction(const DefAmrIndexUint i_input_level,
        const std::set<DefSizet>& surface_process,
        const std::set<std::pair<DefAmrIndexUint, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_reconstruct);
    void ReconstructSurfaceBasedOnExistingVertex(const DefAmrIndexUint i_input_level,
        const std::set<DefSizet>& surface_reconstruct,
        const std::set<std::pair<DefAmrIndexUint, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_remove);
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeometryInfoConnection2D : public GeometryInfo2DInterface, public GeometryConnectionInterface {
 public:
    DefAmrIndexUint k0UxIndex_ = 0, k0UyIndex_ = 0;
    DefAmrIndexUint k0FxIndex_ = 0, k0FyIndex_ = 0;

    // virtual functions for GeometryInfo2DInterface
    void SetIndex() override;
    int InitialGeometry(const DefReal dx, const DefaultGeoShapeType shape_type,
     const DefaultGeoManager& default_geo_managerr) override;
    int UpdateGeometry(const DefaultGeoManager& default_geo_manager) override;
    void FindTrackingNodeBasedOnGeo(const SFBitsetAuxInterface* ptr_sfbitset_aux,
      GridInfoInterface* const ptr_grid_info) override {
        GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(i_geo_, i_level_,
         grid_extend_type_, *ptr_tracking_grid_info_creator_, *ptr_sfbitset_aux, ptr_grid_info);
    }

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel(
        std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max) override;
    DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex1) override;
    void ComputeMidCoordinates(
        const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex1,
        std::vector<DefReal>* const ptr_coordinates) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex] };
    }
    DefSizet GetNumOfGeometryPoints() const final {
        return coordinate_origin_.size();
    }
};
class GeometryInfoConnection2DCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo() override {
        std::shared_ptr<GeometryInfoConnection2D> ptr_temp =
         std::make_shared<GeometryInfoConnection2D>();
        ptr_temp->node_type_ = "Connection2D";
        ptr_temp->geometry_cell_type_ = EGeometryCellType::kPolyLine;
        return ptr_temp;
    };
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeometryInfoConnection3D : public GeometryInfo3DInterface, public GeometryConnectionInterface {
 public:
    DefAmrIndexUint k0UxIndex_ = 0, k0UyIndex_ = 0, k0UzIndex_ = 0;
    DefAmrIndexUint k0FxIndex_ = 0, k0FyIndex_ = 0, k0FzIndex_ = 0;

    // virtual functions for GeometryInfo3DInterface
    void SetIndex() override;
    int InitialGeometry(const DefReal dx,
        const DefaultGeoShapeType shape_type,
        const DefaultGeoManager& default_geo_manager) override;
    int UpdateGeometry(
        const DefaultGeoManager& default_geo_manager) override;
    void FindTrackingNodeBasedOnGeo(const SFBitsetAuxInterface* ptr_sfbitset_aux,
      GridInfoInterface* const ptr_grid_info) override {
        GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(i_geo_, i_level_,
         grid_extend_type_, *ptr_tracking_grid_info_creator_, *ptr_sfbitset_aux, ptr_grid_info);
    }

    // virtual functions for GeometryInfoConnectionInterface
    void InitialCoordinateGivenLevel(
        std::vector<DefReal>* const ptr_coordi_min,
        std::vector<DefReal>* const ptr_coordi_max) override;
    DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex1) override;
    void ComputeMidCoordinates(
        const std::pair<DefAmrIndexUint, DefSizet>& vertex0,
        const std::pair<DefAmrIndexUint, DefSizet>& vertex1,
        std::vector<DefReal>* const ptr_coordinates) override;
    std::vector<DefReal> GetFloodFillOriginArrAsVec() const final {
        return { flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex],
        flood_fill_origin_[kZIndex] };
    }
    DefSizet GetNumOfGeometryPoints() const final {
        return coordinate_origin_.size();
    }
};
class GeometryInfoConnection3DCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo() override {
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
