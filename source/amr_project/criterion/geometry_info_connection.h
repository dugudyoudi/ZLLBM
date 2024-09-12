//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.h
* @author Zhengliang Liu
* @date  2022-9-25
* @brief  define classes to store geometry information with
*         connection relationship of vertices.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_CONNECTION_H_
#define SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_CONNECTION_H_
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
    std::vector<std::pair<DefInt, DefSizet>> vertex_connection;
};
struct GeometryConnectionSurfaceLevel {
    DefInt level_diff;
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
    DefInt level_diff;
    // the first pair<DefInt, DefSizet> is the vertex whose index is larger
    std::map<std::pair<std::pair<DefInt, DefSizet>,
        std::pair<DefInt, DefSizet>>, GeometryConnectionEdge> map_edge_connection;
};
/**
* @struct GeometryVertex
* @brief structure to store vertex formation with connections
*/
struct GeometryVertexConnection : public GeometryVertex {
 public:
    std::map<DefInt, std::set<std::pair<DefInt, DefSizet>>> map_linked_vertices_level;
    ///< indices of vertices at current or different levels linked to this vertex
    std::array<std::pair<DefInt, DefSizet>, 2> parent_vertices;
    std::set<std::pair<DefInt, DefSizet>> child_vertices;
    std::map<DefInt, DefSFBitset> map_bitset_ref;  ///< spacing filing code of vertices
    DefInt highest_grid_level = 0;
    void SetValues(const DefInt level_in, const std::array<std::pair<DefInt, DefSizet>, 2>& parent_in,
        const std::array<DefReal, 3>& coordinate_in) {
        highest_grid_level = level_in;
        parent_vertices = parent_in;
        coordinate = coordinate_in;
    }
};
/**
* @struct GeometryConnectionCoordinateLevel
* @brief structure to store all vertex formation for a given level
*/
struct GeometryConnectionVertexLevel {
    std::vector<std::unique_ptr<GeometryVertexConnection>> vec_vertex_coordinate;
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

    // if true, vertex information other than coordinates,velocity,
    // and force will be stored in vec_int and vec_real for each vertex
    // in GeometryConnectionCoordinate::vertex_info
    bool bool_vertex_info_stored_for_connection_ = false;
    DefInt k0NumEdgeForSurface_ = 3;

    std::vector<std::vector<DefSizet>> connection_relation_;

    std::vector<GeometryConnectionVertexLevel>
        vertex_given_level_{};   ///< vertices at the given level (the ith element)
    std::vector<GeometryConnectionEdgeLevel>
        connection_edge_given_level_{};  ///< edges at the given level (the ith element)
    std::vector<GeometryConnectionSurfaceLevel>
        connection_surface_given_level_{};  ///< surfaces at the given level (ith element)
    // connection_vertex_given_level_ records the vertices at a given level since vertices may exist simultaneously
    // at lower and higher levels, while vertex_given_level_ only store vertices at the lowest levels
    // to reduce memory cost
    std::vector<std::set<std::pair<DefInt, DefSizet>>> connection_vertex_given_level_{};
    ///<  vertices at the current and higher geometry levels exist simultaneously at the given level (the ith element)

    virtual std::unique_ptr<GeometryVertexConnection> GeoConnectionVertexCreator() {
        return std::make_unique<GeometryVertexConnection>();
    }
    virtual std::unique_ptr<GeometryVertexConnection> GeoConnectionVertexCreator(GeometryVertexConnection&) {
        return std::make_unique<GeometryVertexConnection>();
    }

    void InitialCoordinateGivenLevel(const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices,
        std::array<DefReal, 3>* const ptr_coordi_min, std::array<DefReal, 3>* const ptr_coordi_max);
    DefReal ComputeDistanceFromCoordinates(
        const std::pair<DefInt, DefSizet>& vertex0, const std::pair<DefInt, DefSizet>& vertex1);
    void ComputeMidCoordinates(const std::pair<DefInt, DefSizet>& vertex0,
        const std::pair<DefInt, DefSizet>& vertex1, std::array<DefReal, 3>* const ptr_coordinates);

    void SetupConnectionParameters(EGeometryCellType cell_type);
    void InitialConnection(const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices,
        std::array<DefReal, 3>* const ptr_coordi_min, std::array<DefReal, 3>* const ptr_coordi_max);
    void MergeEdgeOnce(const DefInt i_level, const DefInt i_input_level, const DefReal ds_min,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::set<std::pair<std::pair<DefInt, DefSizet>,
        std::pair<DefInt, DefSizet>>>& edge_for_merge,
        std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>*
        const ptr_edhe_remain_for_merge, DefMap<DefInt>* const ptr_sfbitset_ref_removed);
    void BisectEdgeOnce(const DefInt i_level, const DefInt i_input_level, const DefReal ds_max,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::set<std::pair<std::pair<DefInt, DefSizet>,
        std::pair<DefInt, DefSizet>>>& edge_for_bisect,
        std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>*
        const ptr_surface_remain_for_bisect, DefMap<DefInt>* const ptr_sfbitset_ref_added);
    void FindTrackingNodeBasedOnGeo(DefInt i_geo, DefInt i_level,
        const EGridExtendType grid_extend_type, const SFBitsetAuxInterface& sfbitset_aux,
        GridInfoInterface* const ptr_grid_info);

    virtual ~GeometryConnectionInterface() {}

 protected:
    void RemoveVertex(const DefInt i_input_level,
        const std::vector<DefReal>& grid_space,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
        DefMap<DefInt>* const ptr_sfbitset_ref_removed);
    void AddNewLinkage(const DefInt i_input_level,
        const std::pair<DefInt, DefSizet>& vertex_new,
        const std::pair<DefInt, DefSizet>& vertex_origin);
    void FindSurfaceForReconstruction(const DefInt i_input_level,
        const std::set<DefSizet>& surface_process,
        const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_reconstruct);
    void ReconstructSurfaceBasedOnExistingVertex(const DefInt i_input_level,
        const std::set<DefSizet>& surface_reconstruct,
        const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
        std::set<DefSizet>* const ptr_surface_remove);
};
class GeometryInfoConnection : public GeometryInfoInterface, public GeometryConnectionInterface {
 public:
    // virtual functions for GeometryInfo2DInterface
    int InitialGeometry(const DefReal dx) override;
    void FindTrackingNodeBasedOnGeo(const SFBitsetAuxInterface& sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) override {
        GeometryConnectionInterface::FindTrackingNodeBasedOnGeo(i_geo_, i_level_,
         grid_extend_type_, sfbitset_aux, ptr_grid_info);
    }

    explicit GeometryInfoConnection(const DefInt dims) : GeometryInfoInterface(dims) {
        this->node_type_ = "GeometryInfoConnection";
    }
    virtual ~GeometryInfoConnection() {}

 private:
    GeometryInfoConnection();
};
class GeometryInfoConnectionCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo(const DefInt dims) override {
        std::shared_ptr<GeometryInfoConnection> ptr_tmp = std::make_shared<GeometryInfoConnection>(dims);
        return ptr_tmp;
    };
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_CONNECTION_H_
