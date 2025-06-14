//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_connection.cpp
* @author Zhengliang Liu
* @brief functions for geometries with connection relations, e.g. stl files
* @date  2022-10-20
*/
#include <algorithm>
#include <queue>
#include <vector>
#include <string>
#include "./auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_connection.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to initialize status of geometries
* @param[in] dx reference spatial step
* @return  0 successful
*/
void GeometryInfoConnection::InitialGeometry(const DefReal dx) {
    this->SetupConnectionParameters(this->geometry_cell_type_);
    this->GeometryInfoInterface::InitialGeometry(dx);
}
/**
* @brief function to setup geometry type related parameters.
* @param[in]  cell_type    type of geometry cell.
*/
void GeometryConnectionInterface::SetupConnectionParameters(
    EGeometryCellType cell_type) {
    // initialization of cell type related information
    switch (cell_type) {
    case EGeometryCellType::kPolyLine:
        bool_periodic_connection_ = false;
        break;
    case EGeometryCellType::kTriangle:
        bool_periodic_connection_ = true;
        break;
    default:
        LogManager::LogWarning("Undefined cell type.");
        break;
    }
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
* @param[in]   vec_vertices     original vertices of the geometry.
* @param[out]  ptr_coordi_min   minimum coordinates of the geometry.
* @param[out]  ptr_coordi_max   maximum coordinates of the geometry.
*/
void GeometryConnectionInterface::InitialCoordinateGivenLevel(
    const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices,
    std::array<DefReal, 3>* const ptr_coordi_min, std::array<DefReal, 3>* const ptr_coordi_max) {
    vertex_given_level_.push_back({});
    connection_vertex_given_level_.push_back({});
    ptr_coordi_min->at(kXIndex) = vec_vertices.at(0)->coordinate_.at(kXIndex);
    ptr_coordi_min->at(kYIndex) = vec_vertices.at(0)->coordinate_.at(kYIndex);
    ptr_coordi_min->at(kZIndex) = vec_vertices.at(0)->coordinate_.at(kZIndex);
    ptr_coordi_max->at(kXIndex) = vec_vertices.at(0)->coordinate_.at(kXIndex);
    ptr_coordi_max->at(kYIndex) = vec_vertices.at(0)->coordinate_.at(kYIndex);
    ptr_coordi_max->at(kZIndex) = vec_vertices.at(0)->coordinate_.at(kZIndex);
    (vertex_given_level_.at(0).vec_vertex_coordinate_).resize(vec_vertices.size());
    DefSizet i_vertex = 0;
    for (const auto& iter_vertex : vec_vertices) {
        if (ptr_coordi_min->at(kXIndex) > iter_vertex->coordinate_.at(kXIndex)) {
            ptr_coordi_min->at(kXIndex) = iter_vertex->coordinate_.at(kXIndex);
        } else if (ptr_coordi_max->at(kXIndex) < iter_vertex->coordinate_.at(kXIndex)) {
            ptr_coordi_max->at(kXIndex) = iter_vertex->coordinate_.at(kXIndex);
        }
        if (ptr_coordi_min->at(kYIndex) > iter_vertex->coordinate_.at(kYIndex)) {
            ptr_coordi_min->at(kYIndex) = iter_vertex->coordinate_.at(kYIndex);
        } else if (ptr_coordi_max->at(kYIndex) < iter_vertex->coordinate_.at(kYIndex)) {
            ptr_coordi_max->at(kYIndex) = iter_vertex->coordinate_.at(kYIndex);
        }
        if (ptr_coordi_min->at(kZIndex) > iter_vertex->coordinate_.at(kZIndex)) {
            ptr_coordi_min->at(kZIndex) = iter_vertex->coordinate_.at(kZIndex);
        } else if (ptr_coordi_max->at(kZIndex) < iter_vertex->coordinate_.at(kZIndex)) {
            ptr_coordi_max->at(kZIndex) = iter_vertex->coordinate_.at(kZIndex);
        }
        vertex_given_level_.at(0).vec_vertex_coordinate_.at(i_vertex) = GeoConnectionVertexCreator();
        vertex_given_level_.at(0).vec_vertex_coordinate_.at(i_vertex)->coordinate_ = iter_vertex->coordinate_;
        connection_vertex_given_level_.at(0).insert({ 0, i_vertex });
        ++i_vertex;
    }
}
/**
* @brief function to compute distance between two vertices (3D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @return  distance.
*/
DefReal GeometryConnectionInterface::ComputeDistanceFromCoordinates(
    const std::pair<DefInt, DefSizet>& vertex0, const std::pair<DefInt, DefSizet>& vertex1) {
    DefReal x_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kXIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate_.at(vertex1.second)->coordinate_.at(kXIndex));
    DefReal y_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kYIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate_.at(vertex1.second)->coordinate_.at(kYIndex));
    DefReal z_dis = std::fabs(vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kZIndex)
        - vertex_given_level_.at(vertex1.first)
        .vec_vertex_coordinate_.at(vertex1.second)->coordinate_.at(kZIndex));
    return sqrt(x_dis * x_dis + y_dis * y_dis + z_dis * z_dis);
}
/**
* @brief function to compute coordinates of the mid point of two vertices (3D).
* @param[in]  vertex0     indices of the first vertex.
* @param[in]  vertex1     indices of the second vertex.
* @param[out]  ptr_coordinates mid point coordinates.
*/
void GeometryConnectionInterface::ComputeMidCoordinates(const std::pair<DefInt, DefSizet>& vertex0,
    const std::pair<DefInt, DefSizet>& vertex1, std::array<DefReal, 3>* const ptr_coordinates) {
    ptr_coordinates->at(kXIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kXIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate_
        .at(vertex1.second)->coordinate_.at(kXIndex)) / 2.;
    ptr_coordinates->at(kYIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kYIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate_
        .at(vertex1.second)->coordinate_.at(kYIndex)) / 2.;
    ptr_coordinates->at(kZIndex) = (vertex_given_level_.at(vertex0.first)
        .vec_vertex_coordinate_.at(vertex0.second)->coordinate_.at(kZIndex)
        + vertex_given_level_.at(vertex1.first).vec_vertex_coordinate_
        .at(vertex1.second)->coordinate_.at(kZIndex)) / 2.;
}
/**
* @brief function to initialize geometry connections.
* @param[in]   vec_vertices     original vertices of the geometry.
* @param[out]  ptr_coordi_min   minimum coordinates of the geometry.
* @param[out]  ptr_coordi_max   maximum coordinates of the geometry.
*/
void GeometryConnectionInterface::InitialConnection(const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices,
    std::array<DefReal, 3>* const ptr_coordi_min, std::array<DefReal, 3>* const ptr_coordi_max) {
    // step 1: add all geometry points to vertex_given_level_
    InitialCoordinateGivenLevel(vec_vertices, ptr_coordi_min, ptr_coordi_max);
    // add connection relations for edges and surfaces
    GeometryConnectionEdge edge_tmp;
    GeometryConnectionSurface surface_tmp;
    std::pair<DefInt, DefSizet> vertex_index0, vertex_index1;
    DefSizet max_vertex, i_surface = 0, i_vertex;
    connection_edge_given_level_.push_back({});
    connection_edge_given_level_.at(0).level_diff_ = 0;
    connection_surface_given_level_.push_back({});
    connection_surface_given_level_.at(0).level_diff_ = 0;
    connection_vertex_given_level_.push_back({});

    for (const auto& iter_connection : connection_relation_) {
        edge_tmp.set_index_surfaces_.clear();
        edge_tmp.set_index_surfaces_.insert(i_surface);
        max_vertex = iter_connection.size() - 1;
        connection_surface_given_level_.at(0).vec_surface_connection_.emplace_back(surface_tmp);
        for (DefSizet i = 0; i <= max_vertex; ++i) {
            connection_surface_given_level_.at(0).vec_surface_connection_
                .back().vertex_connection_.push_back(
                    std::make_pair(0, iter_connection.at(i)));
            connection_vertex_given_level_.at(0).insert(
                std::make_pair(0, iter_connection.at(i)));
            if (i == max_vertex) {
                if (!bool_periodic_connection_) {
                    break;
                }
                i_vertex = 0;
            } else {
                i_vertex = i + 1;
            }
            if (vertex_given_level_.at(0).vec_vertex_coordinate_
                .at(iter_connection.at(i))->map_linked_vertices_level_.empty()) {
                vertex_given_level_.at(0).vec_vertex_coordinate_
                    .at(iter_connection.at(i))->map_linked_vertices_level_.insert(
                    std::make_pair(0, std::set<std::pair<DefInt, DefSizet>>{}));
            }
            vertex_given_level_.at(0).vec_vertex_coordinate_
                .at(iter_connection.at(i))->map_linked_vertices_level_
                .at(0).insert({ 0, iter_connection.at(i_vertex) });
            if (vertex_given_level_.at(0).vec_vertex_coordinate_
                .at(iter_connection.at(i_vertex))->map_linked_vertices_level_.empty()) {
                vertex_given_level_.at(0).vec_vertex_coordinate_
                    .at(iter_connection.at(i_vertex))->map_linked_vertices_level_.insert(
                    std::make_pair(0, std::set<std::pair<DefInt, DefSizet>>{}));
            }
            vertex_given_level_.at(0).vec_vertex_coordinate_
                .at(iter_connection.at(i_vertex))->map_linked_vertices_level_.at(0).insert(
                { 0, iter_connection.at(i)});
            // edges
            if (iter_connection.at(i) > iter_connection.at(i_vertex)) {
                vertex_index0 = std::make_pair(0, iter_connection.at(i));
                vertex_index1 = std::make_pair(0, iter_connection.at(i_vertex));
            } else {
                vertex_index0 = std::make_pair(0, iter_connection.at(i_vertex));
                vertex_index1 = std::make_pair(0, iter_connection.at(i));
            }
            if (connection_edge_given_level_.at(0).map_edge_connection_
                .find(std::make_pair(vertex_index0, vertex_index1))
                == connection_edge_given_level_.at(0).map_edge_connection_.end()) {
                connection_edge_given_level_.at(0).map_edge_connection_
                    .insert({ std::make_pair(vertex_index0, vertex_index1), edge_tmp });
            } else {
                connection_edge_given_level_.at(0).map_edge_connection_
                    .at(std::make_pair(vertex_index0, vertex_index1))
                    .set_index_surfaces_.insert(i_surface);
            }
        }
        ++i_surface;
    }
}
/**
* @brief function to bisect edges once if they are beyond the threshold.
* @param[in]  i_level    level of geometry.
* @param[in]  i_input_level    level of current connection relations.
* @param[in]  ds_max    upper threshold of the edge length.
* @param[in]  sfbitset_aux class manage space filling curves.
* @param[in]  edge_for_bisect    edges need to be bisected.
* @param[out]  ptr_surface_remain_for_bisect    edge remain to be bisected.
* @param[out]  ptr_sfbitset_ref_added space filling code corresponding to added vertices.            
*/
void GeometryConnectionInterface::BisectEdgeOnce(
    const DefInt i_level, const DefInt i_input_level, const DefReal ds_max,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>& edge_for_bisect,
    std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>*
    const ptr_surface_remain_for_bisect, DefMap<DefInt>* const ptr_sfbitset_ref_added) {
    DefSizet num_surface, num_surface_p1, max_vertex, i_vertex;
    DefReal dis;
    GeometryConnectionEdge edge_tmp;
    GeometryConnectionSurface surface_tmp;
    std::pair<DefInt, DefSizet> vertex_index_tmp, vertex_index_origin;
    std::pair<std::pair<DefInt, DefSizet>,
        std::pair<DefInt, DefSizet>> edge_index_tmp0, edge_index_tmp1;
    // grid space
    std::vector<DefReal> grid_space;
    DefInt i_grid_level = i_level + i_input_level;
    DefReal grid_scale = DefReal(TwoPowerN(i_grid_level));
    DefSFBitset sfbitset_tmp;
    std::vector<DefReal> grid_space_background = sfbitset_aux.GetBackgroundGridSpacing();
    for (const auto iter : grid_space_background) {
        grid_space.emplace_back(iter / grid_scale);
    }
    std::array<DefReal, 3> coordinate;
    std::array<std::pair<DefInt, DefSizet>, 2> parent_vertices_tmp;
    DefInt highest_grid_level = i_grid_level;
    for (const auto& iter_edge : edge_for_bisect) {
        if (connection_edge_given_level_.at(i_input_level).map_edge_connection_.find(iter_edge)
         == connection_edge_given_level_.at(i_input_level).map_edge_connection_.end()) {
            std::string msg = "Can't find edge (";
            msg += std::to_string(iter_edge.first.first) + ", "
                + std::to_string(iter_edge.first.second) + "); ("
                + std::to_string(iter_edge.second.first) + ", "
                + std::to_string(iter_edge.second.second) + ") at level "
                + std::to_string(i_input_level) + " for bisecting.";
            LogManager::LogWarning(msg);
            continue;
        }
        dis = ComputeDistanceFromCoordinates(iter_edge.first, iter_edge.second);
        if (dis > ds_max) {  // edge needs to bisect
            ComputeMidCoordinates(iter_edge.first, iter_edge.second, &coordinate);
            parent_vertices_tmp.at(0) = iter_edge.first;
            parent_vertices_tmp.at(1) = iter_edge.second;
            // indices of the added vertex
            if (iter_edge.first.first > iter_edge.second.first) {
                vertex_index_tmp.first = iter_edge.first.first + 1;
            } else {
                vertex_index_tmp.first = iter_edge.second.first + 1;
            }
            if (static_cast<DefInt>(vertex_given_level_.size()) <= vertex_index_tmp.first) {
                vertex_given_level_.push_back({});
            }
            vertex_given_level_.at(vertex_index_tmp.first)
                .vec_vertex_coordinate_.emplace_back(GeoConnectionVertexCreator());
            vertex_given_level_.at(vertex_index_tmp.first).vec_vertex_coordinate_
                .back()->SetValues(highest_grid_level, parent_vertices_tmp, coordinate);
            // sfbitset of at vertex at i_grid_level
            sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(
                grid_space, {coordinate.at(0), coordinate.at(1), coordinate.at(2)});
            vertex_given_level_.at(vertex_index_tmp.first)
                .vec_vertex_coordinate_.back()->map_bitset_ref_
                .insert({ i_grid_level, sfbitset_tmp });
            if (ptr_sfbitset_ref_added->find(sfbitset_tmp) == ptr_sfbitset_ref_added->end()) {
                ptr_sfbitset_ref_added->insert({ sfbitset_tmp, 1 });
            } else {
                ++ptr_sfbitset_ref_added->at(sfbitset_tmp);
            }
            // child relation
            vertex_index_tmp.second = vertex_given_level_
                .at(vertex_index_tmp.first).vec_vertex_coordinate_.size() - 1;
            vertex_given_level_.at(iter_edge.first.first)
                .vec_vertex_coordinate_.at(iter_edge.first.second)
                ->child_vertices_.insert(vertex_index_tmp);
            vertex_given_level_.at(iter_edge.second.first)
                .vec_vertex_coordinate_.at(iter_edge.second.second)
                ->child_vertices_.insert(vertex_index_tmp);
            if (static_cast<DefInt>(connection_vertex_given_level_.size()) <= i_input_level) {
                connection_vertex_given_level_.push_back({});
            }
            connection_vertex_given_level_.at(i_input_level)
                .insert(vertex_index_tmp);
            // update connection of linked vertices
            vertex_given_level_.at(iter_edge.first.first)
                .vec_vertex_coordinate_.at(iter_edge.first.second)
                ->map_linked_vertices_level_.at(i_input_level).erase(iter_edge.second);
            vertex_given_level_.at(iter_edge.first.first)
                .vec_vertex_coordinate_.at(iter_edge.first.second)
                ->map_linked_vertices_level_.at(i_input_level).insert(vertex_index_tmp);
            vertex_given_level_.at(iter_edge.second.first)
                .vec_vertex_coordinate_.at(iter_edge.second.second)
                ->map_linked_vertices_level_.at(i_input_level).erase(iter_edge.first);
            vertex_given_level_.at(iter_edge.second.first)
                .vec_vertex_coordinate_.at(iter_edge.second.second)
                ->map_linked_vertices_level_.at(i_input_level).insert(vertex_index_tmp);
            // create new surfaces and update connection information
            num_surface = connection_surface_given_level_.
                at(i_input_level).vec_surface_connection_.size();
            for (const auto& iter_surface : connection_edge_given_level_
                .at(i_input_level).map_edge_connection_.at(iter_edge)
                .set_index_surfaces_) {
                num_surface_p1 = num_surface + 1;
                // create new surface
                surface_tmp.parent_surface_ = iter_surface;
                surface_tmp.vertex_connection_ =
                    connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(iter_surface).vertex_connection_;
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.emplace_back(surface_tmp);
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.emplace_back(surface_tmp);
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(iter_surface)
                    .child_surface_.emplace_back(num_surface);
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(iter_surface)
                    .child_surface_.emplace_back(num_surface + 1);

                // update vertex connection
                max_vertex = connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_.at(iter_surface)
                    .vertex_connection_.size() - 1;
                for (DefSizet i = 0; i <= max_vertex; ++i) {
                    // update added vertex of the bisected surface
                    if (connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(iter_surface)
                        .vertex_connection_.at(i) == iter_edge.first) {
                        connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(num_surface)
                            .vertex_connection_.at(i) = vertex_index_tmp;
                    } else if (connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(iter_surface).vertex_connection_.at(i)
                        == iter_edge.second) {
                        connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(num_surface_p1)
                            .vertex_connection_.at(i) = vertex_index_tmp;
                    }

                    // update surface information stored in edges
                    if (i == max_vertex) {
                        if (!bool_periodic_connection_) {
                            break;
                        }
                        i_vertex = 0;
                    } else {
                        i_vertex = i + 1;
                    }
                    // erase parent surface stored in edges
                    if (connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(iter_surface)
                        .vertex_connection_.at(i) >
                        connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(iter_surface)
                        .vertex_connection_.at(i_vertex)) {
                        edge_index_tmp0 = std::make_pair(
                            connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(iter_surface)
                            .vertex_connection_.at(i),
                            connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(iter_surface)
                            .vertex_connection_.at(i_vertex));
                    } else {
                        edge_index_tmp0 = std::make_pair(
                            connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(iter_surface)
                            .vertex_connection_.at(i_vertex),
                            connection_surface_given_level_.at(i_input_level)
                            .vec_surface_connection_.at(iter_surface)
                            .vertex_connection_.at(i));
                    }
                    if (edge_index_tmp0 != iter_edge) {
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.at(edge_index_tmp0)
                            .set_index_surfaces_.erase(iter_surface);
                    }
                    // add bisected surfaces in edges
                    edge_index_tmp1 = edge_index_tmp0;
                    if (connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(iter_surface)
                        .vertex_connection_.at(i) == iter_edge.first) {
                        vertex_index_origin = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(iter_surface).vertex_connection_.at(i_vertex);
                        if (vertex_index_tmp > vertex_index_origin) {
                            edge_index_tmp0 = std::make_pair(
                                vertex_index_tmp, vertex_index_origin);
                        } else {
                            edge_index_tmp0 = std::make_pair(
                                vertex_index_origin, vertex_index_tmp);
                        }
                        AddNewLinkage(i_input_level,
                            vertex_index_tmp, vertex_index_origin);
                    } else if (connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(iter_surface).vertex_connection_.at(i_vertex)
                        == iter_edge.first) {
                        vertex_index_origin = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(iter_surface).vertex_connection_.at(i);
                        if (vertex_index_tmp > vertex_index_origin) {
                            edge_index_tmp0 = std::make_pair(
                                vertex_index_tmp, vertex_index_origin);
                        } else {
                            edge_index_tmp0 = std::make_pair(
                                vertex_index_origin, vertex_index_tmp);
                        }
                        AddNewLinkage(i_input_level,
                            vertex_index_tmp, vertex_index_origin);
                    }
                    if (connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.find(edge_index_tmp0)
                        == connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.end()) {
                        edge_tmp.set_index_surfaces_.clear();
                        edge_tmp.set_index_surfaces_.insert(num_surface);
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.insert(
                                { edge_index_tmp0, edge_tmp });
                        dis = ComputeDistanceFromCoordinates(
                             edge_index_tmp0.first, edge_index_tmp0.second);
                        if (dis > ds_max) {
                            ptr_surface_remain_for_bisect->insert(
                                edge_index_tmp0);
                        }
                    } else {
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.at(edge_index_tmp0)
                            .set_index_surfaces_.insert(num_surface);
                    }
                    if (connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(iter_surface)
                        .vertex_connection_.at(i) == iter_edge.second) {
                        vertex_index_origin = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(iter_surface).vertex_connection_.at(i_vertex);
                        if (vertex_index_tmp > vertex_index_origin) {
                            edge_index_tmp1 = std::make_pair(
                                vertex_index_tmp, vertex_index_origin);
                        } else {
                            edge_index_tmp1 = std::make_pair(
                                vertex_index_origin, vertex_index_tmp);
                        }
                        AddNewLinkage(i_input_level,
                            vertex_index_tmp, vertex_index_origin);
                    } else if (connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(iter_surface).vertex_connection_.at(i_vertex)
                        == iter_edge.second) {
                        vertex_index_origin = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(iter_surface).vertex_connection_.at(i);
                        if (vertex_index_tmp > vertex_index_origin) {
                            edge_index_tmp1 = std::make_pair(
                                vertex_index_tmp, vertex_index_origin);
                        } else {
                            edge_index_tmp1 = std::make_pair(
                                vertex_index_origin, vertex_index_tmp);
                        }
                        AddNewLinkage(i_input_level,
                            vertex_index_tmp, vertex_index_origin);
                    }
                    if (connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.find(edge_index_tmp1)
                        == connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.end()) {
                        edge_tmp.set_index_surfaces_.clear();
                        edge_tmp.set_index_surfaces_.insert(num_surface_p1);
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.insert(
                                { edge_index_tmp1, edge_tmp });
                        dis = ComputeDistanceFromCoordinates(
                            edge_index_tmp1.first, edge_index_tmp1.second);
                        if (dis > ds_max) {
                            ptr_surface_remain_for_bisect->insert(
                                edge_index_tmp1);
                        }
                    } else {
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.at(edge_index_tmp1)
                            .set_index_surfaces_.insert(num_surface_p1);
                    }
                }
                num_surface += 2;
            }
            connection_edge_given_level_.at(i_input_level)
                .map_edge_connection_.erase(iter_edge);
        }
    }
}
/**
* @brief function to merge edge.
* @param[in]  i_level    level of geometry.
* @param[in]  i_input_level    level of current connection relations.
* @param[in]  ds_min    lower threshold of the edge length.
* @param[in]  sfbitset_aux class manage space filling curves.
* @param[in]  edge_for_merge    edges need to be merged.
* @param[out]  ptr_surface_remain_for_bisect    edge remain to be bisected.
* @param[out]  ptr_sfbitset_ref_added space filling code corresponding to added vertices.            
*/
void GeometryConnectionInterface::MergeEdgeOnce(const DefInt i_level,
    const DefInt i_input_level, const DefReal ds_min, const SFBitsetAuxInterface& sfbitset_aux,
    const std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>& edge_for_merge,
    std::set<std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>>*
     const ptr_edge_remain_for_merge, DefMap<DefInt>* const ptr_sfbitset_ref_removed) {
    DefSizet i_vertex, i_vertex0, i_vertex1;
    DefReal dis;
    DefInt base_level = 0;
    DefInt level_vertex0, level_vertex1, level_remove;
    DefSizet vertex_remove;
    std::set<DefSizet> surface_process;
    std::pair<DefInt, DefSizet> vertex_tmp0, vertex_tmp1, vertex_tmp2;
    // (vertex_processed) include (set_vertex_remove) and vertices linked
    // to edges at levels other than i_input_level
    std::set<std::pair<DefInt, DefSizet>> set_vertex_processed,
        set_vertex_remain;
    // find coordinates and edges may need to be removed according
    // to the input
    for (const auto& iter_edge : edge_for_merge) {
        if (connection_edge_given_level_.at(i_input_level).map_edge_connection_.find(iter_edge)
            == connection_edge_given_level_.at(i_input_level).map_edge_connection_.end()) {
            std::string msg = "Can't find edge (" + std::to_string(iter_edge.first.first) + ", "
                + std::to_string(iter_edge.first.second) + "); ("
                + std::to_string(iter_edge.second.first) + ", "
                + std::to_string(iter_edge.second.second) + ") at level "
                + std::to_string(i_input_level) + " for merging.";
            LogManager::LogWarning(msg);
            continue;
        }
        level_vertex0 = iter_edge.first.first;
        i_vertex0 = iter_edge.first.second;
        level_vertex1 = iter_edge.second.first;
        i_vertex1 = iter_edge.second.second;
        dis = ComputeDistanceFromCoordinates(iter_edge.first, iter_edge.second);
        if (dis < ds_min) {  // edge needs to merge
            if (level_vertex0 > level_vertex1) {
                level_remove = level_vertex0;
                vertex_remove = i_vertex0;
            } else if (level_vertex0 < level_vertex1) {
                level_remove = level_vertex1;
                vertex_remove = i_vertex1;
            } else {  // two vertices are at the same layer_level
                if (level_vertex0 == base_level) {
                    continue;
                } else if (i_vertex1 > i_vertex0) {
                    vertex_remove = i_vertex1;
                } else {
                    vertex_remove = i_vertex0;
                }
                level_remove = level_vertex0;
            }
            // coordinates may need to be removed
            vertex_tmp0 = std::make_pair(level_remove, vertex_remove);
            if (set_vertex_processed.find(vertex_tmp0) == set_vertex_processed.end()) {
                set_vertex_processed.insert(vertex_tmp0);
                set_vertex_remain.insert(vertex_tmp0);
                // find surface needs to be reconstructed
                for (const auto& iter_surface : connection_edge_given_level_
                 .at(i_input_level).map_edge_connection_.at(std::make_pair(
                 vertex_tmp0, vertex_given_level_.at(level_remove).vec_vertex_coordinate_
                 .at(vertex_remove)->parent_vertices_.at(0))).set_index_surfaces_) {
                    surface_process.insert(connection_surface_given_level_
                     .at(i_input_level).vec_surface_connection_.at(iter_surface).parent_surface_);
                }
                vertex_tmp1 = vertex_given_level_.at(level_remove)
                 .vec_vertex_coordinate_.at(vertex_remove)->parent_vertices_.at(0);
                vertex_tmp2 = vertex_given_level_.at(level_remove)
                 .vec_vertex_coordinate_.at(vertex_remove)->parent_vertices_.at(1);

                // insert merged edge
                if (vertex_tmp1 >vertex_tmp2) {
                    connection_edge_given_level_.at(i_input_level).map_edge_connection_.
                     insert({ std::make_pair(vertex_tmp1, vertex_tmp2), {} });
                    if (dis / 2. < ds_min) {
                        ptr_edge_remain_for_merge->insert(std::make_pair(vertex_tmp1, vertex_tmp2));
                    }
                } else {
                    connection_edge_given_level_.at(i_input_level).map_edge_connection_
                     .insert({ std::make_pair(vertex_tmp2, vertex_tmp1), {} });
                    if (dis / 2. < ds_min) {
                        ptr_edge_remain_for_merge->insert(std::make_pair(vertex_tmp2, vertex_tmp1));
                    }
                }
            }
        }
    }
    // vertices at current and lower layer levels and their linked edges
    // need to be removed
    set_vertex_processed.clear();
    DefInt iter_count = 0, iter_max = 10;
    std::set<std::pair<DefInt, DefSizet>> set_vertex_remove;

    while (set_vertex_remain.size() > 0 && iter_count < iter_max) {
        ++iter_count;
        std::set<std::pair<DefInt, DefSizet>> set_vertex_remain_tmp
            (set_vertex_remain);
        set_vertex_remain.clear();
        for (const auto& iter_vertex : set_vertex_remain_tmp) {
            if (set_vertex_processed.find(iter_vertex) == set_vertex_processed.end()) {
                set_vertex_processed.insert(iter_vertex);
                // add its child vertices to set_vertex_remain
                for (const auto& iter_child : vertex_given_level_.at(iter_vertex.first)
                    .vec_vertex_coordinate_.at(iter_vertex.second)->child_vertices_) {
                    set_vertex_remain.insert(std::make_pair(iter_child.first, iter_child.second));
                }
                if (vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                    .at(iter_vertex.second)->map_linked_vertices_level_.find(i_input_level)
                    != vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                    .at(iter_vertex.second)->map_linked_vertices_level_.end()) {
                    // delete linkage of the current vertex from others
                    for (const auto& iter_linked_vertex :
                        vertex_given_level_.at(iter_vertex.first)
                        .vec_vertex_coordinate_.at(iter_vertex.second)
                        ->map_linked_vertices_level_.at(i_input_level)) {
                        vertex_given_level_.at(iter_linked_vertex.first)
                            .vec_vertex_coordinate_.at(iter_linked_vertex.second)
                            ->map_linked_vertices_level_.at(i_input_level)
                            .erase(std::make_pair(iter_vertex.first, iter_vertex.second));
                        // remove edges
                        if (iter_vertex > iter_linked_vertex) {
                            connection_edge_given_level_.at(i_input_level).map_edge_connection_
                             .erase(std::make_pair(iter_vertex, iter_linked_vertex));
                        } else {
                            connection_edge_given_level_.at(i_input_level).map_edge_connection_
                             .erase(std::make_pair(iter_linked_vertex, iter_vertex));
                        }
                    }
                    // remove linkage of the current vertex at (i_input_level)
                    vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->map_linked_vertices_level_.erase(i_input_level);
                    if (vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->map_linked_vertices_level_.empty()) {
                        set_vertex_remove.insert(iter_vertex);
                    }
                } else {
                    set_vertex_remove.insert(iter_vertex);
                }
            }
        }
    }
    if (iter_count == iter_max) {
        LogManager::LogWarning("iteration exceeds the maximum in MergeOnceLine");
    }

    // find surface of which vertices need to be removed
    std::set<DefSizet> surface_reconstruct;
    FindSurfaceForReconstruction(i_input_level,
        surface_process, set_vertex_remove, &surface_reconstruct);
    surface_process.clear();

    // reconstruct surface
    std::set<DefSizet> surface_remove;
    ReconstructSurfaceBasedOnExistingVertex(i_input_level, surface_reconstruct,
        set_vertex_remove, &surface_remove);

    // remove surface from (connection_surface_given_level_)
    DefSizet max_vertex, vec_last;
    for (auto iter_surface = surface_remove.rbegin();
        iter_surface != surface_remove.rend(); ++iter_surface) {
        vec_last = connection_surface_given_level_.at(i_input_level).vec_surface_connection_.size() - 1;
        if (*iter_surface < vec_last) {
            // swap the current surface and the last
            std::swap(connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(*iter_surface),
                connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(vec_last));
            // update child and parent relation
            for (const auto iter_c : connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(*iter_surface).child_surface_) {
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(iter_c).parent_surface_
                    = *iter_surface;
            }
            for (auto& iter_c : connection_surface_given_level_
                .at(i_input_level).vec_surface_connection_.at(
                connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(*iter_surface).parent_surface_).child_surface_) {
                if (iter_c == vec_last) {
                    iter_c = *iter_surface;
                    break;
                }
            }
            // update index of surface in edges
            max_vertex = connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(*iter_surface).vertex_connection_.size() - 1;
            for (DefSizet i = 0; i <= max_vertex; ++i) {
                if (i == max_vertex) {
                    if (!bool_periodic_connection_) {
                        break;
                    }
                    i_vertex = 0;
                } else {
                    i_vertex = i + 1;
                }
                if (connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(*iter_surface).vertex_connection_.at(i) >
                    connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(*iter_surface).vertex_connection_.at(i_vertex)) {
                    vertex_tmp0 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(*iter_surface).vertex_connection_.at(i);
                    vertex_tmp1 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(*iter_surface).vertex_connection_.at(i_vertex);
                } else {
                    vertex_tmp1 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(*iter_surface).vertex_connection_.at(i);
                    vertex_tmp0 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(*iter_surface).vertex_connection_.at(i_vertex);
                }
                if (connection_edge_given_level_.at(i_input_level)
                    .map_edge_connection_.find(
                        { vertex_tmp0, vertex_tmp1 })
                    != connection_edge_given_level_.at(i_input_level)
                    .map_edge_connection_.end()) {
                    connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.at({ vertex_tmp0,
                            vertex_tmp1 }).set_index_surfaces_
                        .erase(vec_last);
                    connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.at({ vertex_tmp0,
                            vertex_tmp1 }).set_index_surfaces_
                        .insert(*iter_surface);
                }
            }
        }
        connection_surface_given_level_.at(i_input_level).vec_surface_connection_.pop_back();
    }
    // grid space
    std::vector<DefReal> grid_space;
    DefInt i_grid_level = i_level + i_input_level;
    DefReal grid_scale = DefReal(TwoPowerN(i_grid_level));
    std::vector<DefReal> grid_space_background = sfbitset_aux.GetBackgroundGridSpacing();
    for (const auto iter : grid_space_background) {
        grid_space.emplace_back(iter / grid_scale);
    }
    // remove coordinates from (vertex_given_level_)
    RemoveVertex(i_input_level, grid_space,
        sfbitset_aux, set_vertex_remove, ptr_sfbitset_ref_removed);
}
/**
* @brief function to update the linkage relation of vertex.
* @param[in]  i_input_level    level of current connection relations.
* @param[in]  vertex_new    newly added vertex.
* @param[in]  vertex_origin    vertex already exists.
*/
void GeometryConnectionInterface::AddNewLinkage(const DefInt i_input_level,
    const std::pair<DefInt, DefSizet>& vertex_new,
    const std::pair<DefInt, DefSizet>& vertex_origin) {
    if (vertex_given_level_.at(vertex_new.first)
        .vec_vertex_coordinate_.at(vertex_new.second)
        ->map_linked_vertices_level_.find(i_input_level)
        == vertex_given_level_.at(vertex_new.first)
        .vec_vertex_coordinate_.at(vertex_new.second)
        ->map_linked_vertices_level_.end()) {
        vertex_given_level_.at(vertex_new.first).vec_vertex_coordinate_
            .at(vertex_new.second)->map_linked_vertices_level_.insert({
               i_input_level, {vertex_origin}});
    } else {
        vertex_given_level_.at(vertex_new.first).vec_vertex_coordinate_
            .at(vertex_new.second)->map_linked_vertices_level_.at(i_input_level)
            .insert(vertex_origin);
    }
    vertex_given_level_.at(vertex_origin.first).vec_vertex_coordinate_
        .at(vertex_origin.second)->map_linked_vertices_level_
        .at(i_input_level).insert(vertex_new);
}
/**
* @brief function to remove vertices and their connection relations.
* @param[in]  i_input_level    level of current connection relations.
* @param[in] grid_space   grid spacing
* @param[in]  sfbitset_aux class manage space filling curves.
* @param[in]  set_vertex_remove    vertices need to be removed.
* @param[out] ptr_sfbitset_ref_removed space filling code corresponding
*             to removed vertices.
*/
void GeometryConnectionInterface::RemoveVertex(const DefInt i_input_level,
    const std::vector<DefReal>& grid_space,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
    DefMap<DefInt>* const ptr_sfbitset_ref_removed) {
    DefSizet vec_last;
    std::pair<DefInt, DefSizet> vertex_tmp0, vertex_tmp1;
    std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>> edge_index0, edge_index1;
    DefSFBitset sfbitset_tmp;
    for (auto iter_vertex = set_vertex_remove.rbegin();
        iter_vertex != set_vertex_remove.rend(); ++iter_vertex) {
        // sfbitset of a vertex at i_grid_level
        std::array<DefReal, 3>& coordi(vertex_given_level_.at(iter_vertex->first)
            .vec_vertex_coordinate_.at(iter_vertex->second)->coordinate_);
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(
            grid_space, {coordi[0], coordi[1], coordi[2]});
        if (ptr_sfbitset_ref_removed->find(sfbitset_tmp) == ptr_sfbitset_ref_removed->end()) {
            ptr_sfbitset_ref_removed->insert({ sfbitset_tmp, 1 });
        } else {
            ++ptr_sfbitset_ref_removed->at(sfbitset_tmp);
        }
        // delete child relation for the vertex to be removed
        vertex_tmp0 = vertex_given_level_.at(iter_vertex->first)
            .vec_vertex_coordinate_.at(iter_vertex->second)->parent_vertices_.at(0);
        vertex_given_level_.at(vertex_tmp0.first)
            .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.erase(*iter_vertex);
        vertex_tmp0 = vertex_given_level_.at(iter_vertex->first)
            .vec_vertex_coordinate_.at(iter_vertex->second)->parent_vertices_.at(1);
        vertex_given_level_.at(vertex_tmp0.first)
            .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.erase(*iter_vertex);

        vec_last = vertex_given_level_.at(iter_vertex->first)
            .vec_vertex_coordinate_.size() - 1;
        if (iter_vertex->second < vec_last) {
            // update child relation for the last vertex
            vertex_tmp1 = std::make_pair(iter_vertex->first, vec_last);
            vertex_tmp0 = vertex_given_level_.at(vertex_tmp1.first)
                .vec_vertex_coordinate_.at(vertex_tmp1.second)->parent_vertices_.at(0);
            vertex_given_level_.at(vertex_tmp0.first)
                .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.erase(vertex_tmp1);
            vertex_given_level_.at(vertex_tmp0.first)
                .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.insert(*iter_vertex);
            vertex_tmp0 = vertex_given_level_.at(vertex_tmp1.first)
                .vec_vertex_coordinate_.at(vertex_tmp1.second)->parent_vertices_.at(1);
            vertex_given_level_.at(vertex_tmp0.first)
                .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.erase(vertex_tmp1);
            vertex_given_level_.at(vertex_tmp0.first)
                .vec_vertex_coordinate_.at(vertex_tmp0.second)->child_vertices_.insert(*iter_vertex);
            // swap vertices
            std::swap(vertex_given_level_.at(iter_vertex->first)
                .vec_vertex_coordinate_.at(iter_vertex->second),
                vertex_given_level_.at(iter_vertex->first)
                .vec_vertex_coordinate_.at(vec_last));
            if (connection_vertex_given_level_.at(i_input_level)
                .find({iter_vertex->first, vec_last})
                != connection_vertex_given_level_.at(i_input_level).end()) {
                connection_vertex_given_level_.at(i_input_level)
                    .erase({ iter_vertex->first, vec_last });
            }
            // update parent relation
            for (auto& iter_child : vertex_given_level_.at(iter_vertex->first)
                .vec_vertex_coordinate_.at(iter_vertex->second)->child_vertices_) {
                if (vertex_given_level_.at(iter_child.first).vec_vertex_coordinate_
                    .at(iter_child.second)->parent_vertices_.at(0).first == iter_vertex->first
                    && vertex_given_level_.at(iter_child.first).vec_vertex_coordinate_
                    .at(iter_child.second)->parent_vertices_.at(0).second == vec_last) {
                    vertex_given_level_.at(iter_child.first).vec_vertex_coordinate_
                        .at(iter_child.second)->parent_vertices_.at(0).second = iter_vertex->second;
                } else {
                    vertex_given_level_.at(iter_child.first).vec_vertex_coordinate_
                        .at(iter_child.second)->parent_vertices_.at(1).second = iter_vertex->second;
                }
            }
            DefInt level_diff;
            for (const auto& iter_link_level : vertex_given_level_
                .at(iter_vertex->first).vec_vertex_coordinate_
                .at(iter_vertex->second)->map_linked_vertices_level_) {
                // at all the levels connected to the vertex
                level_diff = iter_link_level.first;
                for (const auto& iter_link : iter_link_level.second) {
                    vertex_tmp0 = { iter_vertex->first, vec_last };
                    // update linkage
                    vertex_given_level_.at(iter_link.first)
                        .vec_vertex_coordinate_.at(iter_link.second)
                        ->map_linked_vertices_level_.at(level_diff).erase(vertex_tmp0);
                    vertex_given_level_.at(iter_link.first)
                        .vec_vertex_coordinate_.at(iter_link.second)
                        ->map_linked_vertices_level_.at(level_diff).insert(*iter_vertex);
                    // update edge and surface for the swapped node
                    if (vertex_tmp0 > iter_link) {
                        edge_index0 = { vertex_tmp0,  iter_link };
                    } else {
                        edge_index0 = { iter_link, vertex_tmp0};
                    }
                    for (const auto& iter_surface : connection_edge_given_level_
                        .at(level_diff).map_edge_connection_.at(edge_index0).set_index_surfaces_) {
                        for (auto& iter_sf_vertex : connection_surface_given_level_.at(level_diff)
                            .vec_surface_connection_.at(iter_surface).vertex_connection_) {
                            if (iter_sf_vertex == vertex_tmp0) {
                                iter_sf_vertex = *iter_vertex;
                            }
                        }
                    }
                    if (*iter_vertex > iter_link) {
                        edge_index1 = { *iter_vertex,  iter_link };
                    } else {
                        edge_index1 = { iter_link , *iter_vertex};
                    }
                    if (connection_edge_given_level_.at(level_diff)
                        .map_edge_connection_.find(edge_index1)
                        == connection_edge_given_level_.at(level_diff)
                        .map_edge_connection_.end()) {
                        auto edge_move = connection_edge_given_level_
                            .at(level_diff).map_edge_connection_.extract(edge_index0);
                        edge_move.key() = edge_index1;
                        connection_edge_given_level_.at(level_diff)
                            .map_edge_connection_.insert(std::move(edge_move));
                    } else {
                        connection_edge_given_level_.at(level_diff).map_edge_connection_
                            .at(edge_index1).set_index_surfaces_.swap(
                            connection_edge_given_level_.at(level_diff)
                            .map_edge_connection_.at(edge_index0).set_index_surfaces_);
                        connection_edge_given_level_.at(level_diff)
                            .map_edge_connection_.erase(edge_index0);
                    }
                }
            }
        }
        vertex_given_level_.at(iter_vertex->first).vec_vertex_coordinate_.pop_back();
    }
}
/**
* @brief function to find surface need to be reconstructed.
* @param[in]  i_input_level    level of current connection relations.
* @param[in]  surface_process    input surfaces.
* @param[in]  set_vertex_remove    vertices need to be removed.
* @param[out] ptr_surface_reconstruct surfaces need to be reconstructed.
*/
void GeometryConnectionInterface::FindSurfaceForReconstruction(
    const DefInt i_input_level, const std::set<DefSizet>& surface_process,
    const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
    std::set<DefSizet>* const ptr_surface_reconstruct) {
    bool bool_all_exist;
    DefSizet i_surface;
    for (const auto& iter_surface : surface_process) {
        bool_all_exist = false;
        i_surface = iter_surface;
        while (!bool_all_exist) {
            bool_all_exist = true;
            for (auto iter : connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection_.at(i_surface).vertex_connection_) {
                if (set_vertex_remove.find(iter) != set_vertex_remove.end()) {
                    bool_all_exist = false;
                    i_surface = connection_surface_given_level_.at(i_input_level)
                        .vec_surface_connection_.at(i_surface).parent_surface_;
                    break;
                }
            }
        }
        ptr_surface_reconstruct->insert(i_surface);
    }
}
/**
* @brief function to reconstruct surface based on existing vertices.
* @param[in]  i_input_level    level of current connection relations.
* @param[in]  surface_reconstruct    surfaces need to be reconstructed.
* @param[in]  set_vertex_remove    vertices need to be removed.
* @param[out] ptr_surface_remove surfaces need to be removed.
*/
void GeometryConnectionInterface::ReconstructSurfaceBasedOnExistingVertex(
    const DefInt i_input_level, const std::set<DefSizet>& surface_reconstruct,
    const std::set<std::pair<DefInt, DefSizet>>& set_vertex_remove,
    std::set<DefSizet>* const ptr_surface_remove) {
    DefSizet i_surface, i_surface_out;
    GeometryConnectionEdge edge_connection_tmp;
    std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>>
        edge_key_tmp;
    DefSizet i_vertex, max_vertex;
    std::pair<DefInt, DefSizet> vertex_tmp0, vertex_tmp1;
    for (auto iter_surface = surface_reconstruct.rbegin();
        iter_surface != surface_reconstruct.rend(); ++iter_surface) {
        std::queue<DefSizet> surface_tmp, surface_remain;
        std::map<std::pair<std::pair<DefInt, DefSizet>,
            std::pair<DefInt, DefSizet>>, std::pair<DefInt, DefSizet>>
            edge_of_midpoint;
        surface_tmp.push(*iter_surface);

        // find vertices do not need to remove in child surfaces
        while (!surface_tmp.empty()) {
            i_surface = surface_tmp.front();
            surface_tmp.pop();
            surface_remain.push(i_surface);
            if (connection_surface_given_level_
                .at(i_input_level).vec_surface_connection_.at(i_surface)
                .child_surface_.empty()) {
                // delete surface from edges (assuming the vertex are connected
                // sequentially)
                max_vertex = connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_.at(i_surface)
                    .vertex_connection_.size() - 1;
                for (DefSizet i = 0; i <= max_vertex; ++i) {
                    if (i == max_vertex) {
                        if (!bool_periodic_connection_) {
                            break;
                        }
                        i_vertex = 0;
                    } else {
                        i_vertex = i + 1;
                    }
                    if (connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i) >
                        connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i_vertex)) {
                        vertex_tmp0 = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(i_surface).vertex_connection_.at(i);
                        vertex_tmp1 = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(i_surface).vertex_connection_.at(i_vertex);
                    } else {
                        vertex_tmp1 = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(i_surface).vertex_connection_.at(i);
                        vertex_tmp0 = connection_surface_given_level_
                            .at(i_input_level).vec_surface_connection_
                            .at(i_surface).vertex_connection_.at(i_vertex);
                    }
                    if (connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.find(
                            { vertex_tmp0, vertex_tmp1 })
                        != connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.end()) {
                        connection_edge_given_level_.at(i_input_level)
                            .map_edge_connection_.at({ vertex_tmp0,
                                vertex_tmp1 }).set_index_surfaces_
                            .erase(i_surface);
                    }
                }
            } else {
                // push child surfaces
                while (!connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_.at(i_surface)
                    .child_surface_.empty()) {
                    surface_tmp.push(connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_.at(i_surface)
                        .child_surface_.back());
                    connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_.at(i_surface)
                        .child_surface_.pop_back();
                }
            }
            // find child vertices do not need to be removed
            for (const auto& iter_vertex : connection_surface_given_level_
                .at(i_input_level).vec_surface_connection_.at(i_surface).vertex_connection_) {
                if (set_vertex_remove.find(iter_vertex)
                    == set_vertex_remove.end()) {
                    // if parent vertices exist, i.e. parent_vertices[0]
                    // != parent_vertices[1]
                    if (vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->parent_vertices_.at(0) >
                        vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->parent_vertices_.at(1)) {
                        edge_of_midpoint.insert({ std::make_pair(
                            vertex_given_level_.at(iter_vertex.first)
                            .vec_vertex_coordinate_.at(iter_vertex.second)->parent_vertices_.at(0),
                            vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                            .at(iter_vertex.second)->parent_vertices_.at(1)), iter_vertex });
                    } else if (vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->parent_vertices_.at(0)
                        < vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                        .at(iter_vertex.second)->parent_vertices_.at(1)) {
                        edge_of_midpoint.insert({ std::make_pair(
                            vertex_given_level_.at(iter_vertex.first)
                            .vec_vertex_coordinate_.at(iter_vertex.second)->parent_vertices_.at(1),
                            vertex_given_level_.at(iter_vertex.first).vec_vertex_coordinate_
                            .at(iter_vertex.second)->parent_vertices_.at(0)), iter_vertex });
                    }
                }
            }
        }
        // bisect one edge of surfaces for all exist vertices
        surface_tmp.push(*iter_surface);
        bool bool_child;
        DefSizet count_vertex0, count_vertex1;
        std::pair<std::pair<DefInt, DefSizet>, std::pair<DefInt, DefSizet>> edge_key_min;
        while (!surface_tmp.empty()) {
            i_surface = surface_tmp.front();
            surface_tmp.pop();
            edge_key_min = { { ~0, ~0 }, { ~0, ~0 } };
            edge_connection_tmp.set_index_surfaces_ = { i_surface };
            bool_child = false;
            max_vertex = connection_surface_given_level_
                .at(i_input_level).vec_surface_connection_.at(i_surface)
                .vertex_connection_.size() - 1;
            for (DefSizet i = 0; i <= max_vertex; ++i) {
                if (i == max_vertex) {
                    if (!bool_periodic_connection_) {
                        break;
                    }
                    i_vertex = 0;
                } else {
                    i_vertex = i + 1;
                }
                if (connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(i_surface).vertex_connection_.at(i) >
                    connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(i_surface).vertex_connection_.at(i_vertex)) {
                    edge_key_tmp = { connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i),
                        connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i_vertex) };
                } else {
                    edge_key_tmp = { connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i_vertex),
                        connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_
                        .at(i_surface).vertex_connection_.at(i) };
                }
                if (edge_of_midpoint.find(edge_key_tmp) != edge_of_midpoint.end()
                    && edge_key_tmp < edge_key_min) {
                    edge_key_min = edge_key_tmp;
                    count_vertex0 = i;
                    count_vertex1 = i_vertex;
                    bool_child = true;
                }
            }
            if (bool_child) {
                // bisect the surface
                surface_remain.pop();
                i_surface_out = surface_remain.front();
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface)
                    .child_surface_.emplace_back(i_surface_out);
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .parent_surface_ = i_surface;
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .child_surface_.clear();
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .vertex_connection_ = connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(i_surface).vertex_connection_;
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .vertex_connection_.at(count_vertex0)
                    = edge_of_midpoint.at(edge_key_min);
                surface_tmp.push(i_surface_out);

                surface_remain.pop();
                i_surface_out = surface_remain.front();
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface)
                    .child_surface_.emplace_back(i_surface_out);
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .parent_surface_ = i_surface;
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .child_surface_.clear();
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .vertex_connection_ = connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_
                    .at(i_surface).vertex_connection_;
                connection_surface_given_level_.at(i_input_level)
                    .vec_surface_connection_.at(i_surface_out)
                    .vertex_connection_.at(count_vertex1)
                    = edge_of_midpoint.at(edge_key_min);
                surface_tmp.push(i_surface_out);

                edge_of_midpoint.erase(edge_key_min);
            } else {  // surface do not need to reconstruct anymore
                max_vertex = connection_surface_given_level_
                    .at(i_input_level).vec_surface_connection_.at(i_surface)
                    .vertex_connection_.size() - 1;
                for (DefSizet i = 0; i <= max_vertex; ++i) {
                    if (i == max_vertex) {
                        if (!bool_periodic_connection_) {
                            break;
                        }
                        i_vertex = 0;
                    } else {
                        i_vertex = i + 1;
                    }
                    vertex_tmp0 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_.at(i_surface)
                        .vertex_connection_.at(i);
                    vertex_tmp1 = connection_surface_given_level_
                        .at(i_input_level).vec_surface_connection_.at(i_surface)
                        .vertex_connection_.at(i_vertex);

                    // update vertex linkage connection
                    if (vertex_given_level_.at(vertex_tmp0.first)
                        .vec_vertex_coordinate_.at(vertex_tmp0.second)
                        ->map_linked_vertices_level_.find(i_input_level)
                        == vertex_given_level_.at(vertex_tmp0.first)
                        .vec_vertex_coordinate_.at(vertex_tmp0.second)
                        ->map_linked_vertices_level_.end()) {
                        vertex_given_level_.at(vertex_tmp0.first)
                            .vec_vertex_coordinate_.at(vertex_tmp0.second)
                            ->map_linked_vertices_level_.insert({ i_input_level,
                            {vertex_tmp1} });
                    } else {
                        vertex_given_level_.at(vertex_tmp0.first)
                            .vec_vertex_coordinate_.at(vertex_tmp0.second)
                            ->map_linked_vertices_level_.at(i_input_level)
                            .insert(vertex_tmp1);
                    }
                    if (vertex_given_level_.at(vertex_tmp1.first)
                        .vec_vertex_coordinate_.at(vertex_tmp1.second)
                        ->map_linked_vertices_level_.find(i_input_level)
                        == vertex_given_level_.at(vertex_tmp1.first)
                        .vec_vertex_coordinate_.at(vertex_tmp1.second)
                        ->map_linked_vertices_level_.end()) {
                        vertex_given_level_.at(vertex_tmp1.first)
                            .vec_vertex_coordinate_.at(vertex_tmp1.second)
                            ->map_linked_vertices_level_.insert({ i_input_level,
                                {vertex_tmp0} });
                    } else {
                        vertex_given_level_.at(vertex_tmp1.first)
                            .vec_vertex_coordinate_.at(vertex_tmp1.second)
                            ->map_linked_vertices_level_.at(i_input_level)
                            .insert(vertex_tmp0);
                    }
                    // update edge
                    if (vertex_tmp0 > vertex_tmp1) {
                        edge_key_tmp = std::make_pair(vertex_tmp0, vertex_tmp1);
                    } else {
                        edge_key_tmp = std::make_pair(vertex_tmp1, vertex_tmp0);
                    }
                    if (connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.find(edge_key_tmp)
                        == connection_edge_given_level_.at(i_input_level)
                        .map_edge_connection_.end()) {
                        connection_edge_given_level_.at(i_input_level).map_edge_connection_
                            .insert(std::make_pair(edge_key_tmp, edge_connection_tmp));
                    } else {
                        connection_edge_given_level_.at(i_input_level).map_edge_connection_
                            .at(edge_key_tmp).set_index_surfaces_.insert(i_surface);
                    }
                }
            }
        }

        // surface need to be removed
        surface_remain.pop();
        while (!surface_remain.empty()) {
            ptr_surface_remove->insert(surface_remain.front());
            surface_remain.pop();
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
