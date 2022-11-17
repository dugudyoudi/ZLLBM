//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_stl.cpp
* @author Zhengliang Liu
* @brief functions for geometries reprensented by STL points
* @date  2022-5-25
* @note .
*/
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
#include <algorithm>
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "criterion/geometry_info_connection.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
namespace criterion {
/**
* @brief   function to set indices and size of vectors in GeometryCoordinate to
*          store real and int variable
* @note
*/
void GeometryInfoConnection2D::SetIndex() {
    if (bool_vec_velocities_) {
        k0UxIndex_ = 0; k0UyIndex_ = k0UxIndex_ + 1;
        k0NumRealForEachPoint_ = k0UyIndex_;
    }
    if (bool_vec_forces_) {
        k0FxIndex_ = k0UyIndex_ + 1; k0FyIndex_ = k0FxIndex_ + 1;
        k0NumRealForEachPoint_ = k0FyIndex_;
    }

    geo_point_info_instance_.vec_real.resize(k0NumRealForEachPoint_);
}

/**
* @brief   function to initialize status of geometries
* @return  0 successful, 1 geometry is undefined
* @note
*/
int GeometryInfoConnection2D::InitialGeometry(
    std::shared_ptr<GeometryInfo2DInterface> ptr_geo) {
    DefUint dims = 2;
    SetIndex();
    SetupConnectionParameters(geometry_cell_type_);
    switch (k0DefaultGeoShapeType_) {
    case DefaultGeoShapeType::kCircle:
        DefaultGeoShapeManager::GetInstance()
            ->geo_circle_initial(dims, ptr_geo);
        return 0;
    default:
        return 1;
    }
}
int GeometryInfoConnection2D::UpdateGeometry(
    std::shared_ptr<GeometryInfo2DInterface> ptr_geo) {
    return 0;
}
/**
* @brief   function to copy coordinates from vector of original coordinates
*           to that used for identifying connection relations.
*/
void GeometryInfoConnection2D::InitialCoordinateGivenLevel() {
    coordinate_given_level_.push_back({});
    connection_index_given_level_.push_back({});
    GeometryConnectionCoordinate2D coordi_temp;
    coordi_temp.vec_index_edges = {};
    coordi_temp.point_info = connection_vertex_info_instance_;
    DefSizet i_point = 0;
    for (const auto& iter_point : coordinate_origin_) {
        coordi_temp.coordinate = iter_point.coordinate;
        coordinate_given_level_.at(0).vec_vertex_cooridinate
            .push_back(coordi_temp);
        connection_index_given_level_.at(0)
            .insert(std::make_pair(0, i_point));
        ++i_point;
    }
}
void GeometryInfoConnection2D::InitialConnectionAddEdge(
    const DefSizet index_small, const DefSizet index_large,
    const DefSizet i_surface) {
    // check if the edge exists
    bool bool_exist = false;
    if (coordinate_given_level_.at(0).vec_vertex_cooridinate
        .at(index_small).vec_index_edges.size() > 0) {
        for (const auto& iter_edge :
            coordinate_given_level_.at(0).vec_vertex_cooridinate
            .at(index_small).vec_index_edges) {
            if (connection_edge_given_level_.at(0)
                .vec_edge_connection.at(iter_edge.second)
                .vertex_connection[1].vertex_index == index_large) {
                bool_exist = true;
                connection_edge_given_level_.at(0)
                    .vec_edge_connection.at(iter_edge.second)
                    .vec_index_surfaces.push_back(i_surface);
                connection_surface_given_level_.at(0)
                    .vec_surface_connection.at(i_surface)
                    .edge_connection.push_back(iter_edge.second);
                break;
            }
        }
    }

    GeometryConnectionEdge connection_edge_temp;
    GeometryCoordinateIndex vertex_index_temp1, vertex_index_temp2;
    vertex_index_temp1.layer_level = 0;
    vertex_index_temp2.layer_level = 0;
    if (!bool_exist) {
        vertex_index_temp1.vertex_index = index_small;
        vertex_index_temp2.vertex_index = index_large;
        connection_edge_temp.vertex_connection = {
        vertex_index_temp1, vertex_index_temp2 };
        DefSizet num_edge = connection_edge_given_level_
            .at(0).vec_edge_connection.size();
        std::pair<DefSizet, DefSizet> edge_index = std::make_pair(0, num_edge);
        connection_edge_temp.vec_index_surfaces = { i_surface };
        connection_edge_given_level_.at(0).vec_edge_connection
            .push_back(connection_edge_temp);
        coordinate_given_level_.at(0).vec_vertex_cooridinate
            .at(index_small).vec_index_edges.push_back(edge_index);
        coordinate_given_level_.at(0).vec_vertex_cooridinate
            .at(index_large).vec_index_edges.push_back(edge_index);
        connection_surface_given_level_.at(0).vec_surface_connection
            .at(i_surface).edge_connection.push_back(num_edge);
    }
}
void GeometryInfoConnection2D::BisectinOnceLine(
    const DefSizet i_input_level, const DefReal ds_max,
    const std::vector<DefSizet>& vec_edge_for_bisect,
    std::vector<DefSizet>* ptr_surface_remain_for_bisect) {
    DefSizet i_vertex, i_vertex0, i_vertex1, i_layer_level;
    DefSizet num_edge, num_surface;
    DefReal dis;
    GeometryConnectionCoordinate2D coordi_connection_temp;
    GeometryConnectionEdge edge_temp;
    GeometryConnectionSurface surface_temp;
    std::array<DefReal, 2> coordi0, coordi1;
    DefSizet level_vertex0, level_vertex1;
    //if (coordinate_given_level_.size() < i_output_level + 1) {
    //    // add connection information at (i_output_level) level
    //    coordinate_given_level_.push_back({});
    //    coordinate_given_level_.back().i_layer_level
    //        = i_level_ + i_output_level;

    //    connection_edge_given_level_.push_back({});
    //    connection_edge_given_level_.back().i_layer_level
    //        = i_level_ + i_output_level;

    //    connection_surface_given_level_.push_back({});
    //    connection_surface_given_level_.back().i_layer_level
    //        = i_level_ + i_output_level;

    //    connection_index_given_level_.push_back({});
    //}

    for (const auto& iter_edge : vec_edge_for_bisect) {
        level_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).layer_level;
        i_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).vertex_index;
        level_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).layer_level;
        i_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).vertex_index;
        coordi0 = coordinate_given_level_.at(level_vertex0)
            .vec_vertex_cooridinate.at(i_vertex0).coordinate;
        coordi1 = coordinate_given_level_.at(level_vertex1)
            .vec_vertex_cooridinate.at(i_vertex1).coordinate;
        dis = ComputeDistanceFromCoordinates(coordi0, coordi1);
        if (dis > ds_max) {  // edge needs to bisect
            num_edge = connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.size();
            num_surface = connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection.size();
            // add coordinates
            coordi_connection_temp.coordinate.at(kXIndex)
                = (coordi0.at(kXIndex) + coordi1.at(kXIndex)) / 2.;
            coordi_connection_temp.coordinate.at(kYIndex)
                = (coordi0.at(kYIndex) + coordi1.at(kYIndex)) / 2.;
            coordi_connection_temp.parent_vertices.at(0)
                = connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.at(iter_edge).vertex_connection.at(0);
            coordi_connection_temp.parent_vertices.at(1)
                = connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.at(iter_edge).vertex_connection.at(1);
            coordi_connection_temp.vec_index_edges = {
                std::make_pair(i_input_level, iter_edge),
                std::make_pair(i_input_level, num_edge) };
            if (level_vertex0 > level_vertex1) {
                i_layer_level = level_vertex0 + 1;
            } else {
                i_layer_level = level_vertex1 + 1;
            }
            if (coordinate_given_level_.size() > i_layer_level) {
                coordinate_given_level_.at(i_layer_level)
                    .vec_vertex_cooridinate.push_back(coordi_connection_temp);
            } else {
                coordinate_given_level_.push_back({});
                coordinate_given_level_.at(i_layer_level)
                    .vec_vertex_cooridinate = { coordi_connection_temp };
            }
            i_vertex = coordinate_given_level_.at(i_layer_level)
                .vec_vertex_cooridinate.size() - 1;
            coordinate_given_level_.at(level_vertex0).vec_vertex_cooridinate
                .at(i_vertex0).vec_child_vertices
                .push_back(std::make_pair(i_layer_level, i_vertex));
            coordinate_given_level_.at(level_vertex1).vec_vertex_cooridinate
                .at(i_vertex1).vec_child_vertices
                .push_back(std::make_pair(i_layer_level, i_vertex));
            // update connection relation
            for (auto& iter_edge_index : coordinate_given_level_
                .at(level_vertex1).vec_vertex_cooridinate.at(i_vertex1)
                .vec_index_edges) {
                if (iter_edge_index.first == i_input_level
                    && iter_edge_index.second == iter_edge) {
                    iter_edge_index.second = num_edge;
                    break;
                }
            }
            edge_temp.vertex_connection.at(0) =
                connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.at(iter_edge).vertex_connection.at(1);
            connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.at(iter_edge)
                .vertex_connection.at(1).layer_level = i_layer_level;
            connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.at(iter_edge)
                .vertex_connection.at(1).vertex_index = i_vertex;
            edge_temp.vertex_connection.at(1) = connection_edge_given_level_
                .at(i_input_level).vec_edge_connection.at(iter_edge)
                .vertex_connection.at(1);
            edge_temp.vec_index_surfaces = { num_surface };
            connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.push_back(edge_temp);
            surface_temp.edge_connection = { num_edge };
            connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection.push_back(surface_temp);
            connection_index_given_level_.at(i_input_level).insert(
                std::make_pair(i_layer_level, i_vertex));
            if (dis / 2. > ds_max) {
                ptr_surface_remain_for_bisect->push_back(iter_edge);
                ptr_surface_remain_for_bisect->push_back(num_edge);
            }
        }
    }
}
void GeometryInfoConnection2D::MergeOnceLine(
    const DefSizet i_input_level, const DefReal ds_min,
    const std::vector<DefSizet>& vec_edge_for_bisect,
    std::vector<DefSizet>* ptr_surface_remain_for_merge) {
    DefSizet i_vertex, i_vertex0, i_vertex1, i_layer_level;
    DefSizet num_edge, num_surface, edge_index_min;
    DefReal dis;
    GeometryConnectionCoordinate2D coordi_connection_temp;
    GeometryConnectionEdge edge_temp;
    GeometryConnectionSurface surface_temp;
    std::array<DefReal, 2> coordi0, coordi1;
    DefSizet level_vertex0, level_vertex1, level_remove, vertex_remove;
    std::set<DefSizet>  edge_remove, edge_remain, edge_for_merge, surface_move;
    std::pair<DefSizet, DefSizet> point_index_temp;
    // (point_porcessed) include (point_remove) and points linked to edges at
    // levels other than i_input_level
    std::set<std::pair<DefSizet, DefSizet>> point_processed,
        point_remove, point_remain;
    for (const auto& iter_edge : vec_edge_for_bisect) {
        level_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).layer_level;
        i_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).vertex_index;
        level_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).layer_level;
        i_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).vertex_index;
        coordi0 = coordinate_given_level_.at(level_vertex0)
            .vec_vertex_cooridinate.at(i_vertex0).coordinate;
        coordi1 = coordinate_given_level_.at(level_vertex1)
            .vec_vertex_cooridinate.at(i_vertex1).coordinate;
        dis = ComputeDistanceFromCoordinates(coordi0, coordi1);
        if (dis < ds_min) {  // edge needs to merge
            num_edge = connection_edge_given_level_.at(i_input_level)
                .vec_edge_connection.size();
            num_surface = connection_surface_given_level_.at(i_input_level)
                .vec_surface_connection.size();
            if (level_vertex0 > level_vertex1) {
                level_remove = level_vertex0;
                vertex_remove = i_vertex0;
            } else if (level_vertex0 < level_vertex1) {
                level_remove = level_vertex1;
                vertex_remove = i_vertex1;
            } else {  // two vertices are at the same layer_level
                continue;
            }
            // coordinates may need to be removed
            point_index_temp = std::make_pair(level_remove, vertex_remove);
            if (point_processed.find(point_index_temp)
                == point_processed.end()) {
                point_processed.insert(point_index_temp);
                edge_index_min = num_edge + 1;
                std::vector<std::pair<DefSizet, DefSizet>> vec_edge_others{};
                for (const auto& iter_edge_index : coordinate_given_level_
                    .at(level_remove).vec_vertex_cooridinate.at(vertex_remove)
                    .vec_index_edges) {
                    // only take actions for edges at (i_input_level)
                    if (iter_edge_index.first == i_input_level) {
                        if (edge_remove.find(iter_edge_index.second)
                            == edge_remove.end()) {
                            edge_remove.insert(iter_edge_index.second);
                        }
                        if (edge_index_min > iter_edge_index.second) {
                            edge_index_min = iter_edge_index.second;
                        }
                        // find childrens whose (layer_level) is higher than
                        // the (level_remove), which also need to be removed
                        if (connection_edge_given_level_.at(i_input_level)
                            .vec_edge_connection.at(iter_edge_index.second)
                            .vertex_connection.at(0).layer_level
                             > level_remove) {
                            point_remain.insert(std::make_pair(
                                connection_edge_given_level_.at(i_input_level)
                                .vec_edge_connection.at(iter_edge_index.second)
                                .vertex_connection.at(0).layer_level,
                                connection_edge_given_level_.at(i_input_level)
                                .vec_edge_connection.at(iter_edge_index.second)
                                .vertex_connection.at(0).vertex_index));
                        } else if (connection_edge_given_level_
                            .at(i_input_level).vec_edge_connection
                            .at(iter_edge_index.second).vertex_connection.at(1)
                            .layer_level > level_remove) {
                            point_remain.insert(std::make_pair(
                                connection_edge_given_level_.at(i_input_level)
                                .vec_edge_connection.at(iter_edge_index.second)
                                .vertex_connection.at(1).layer_level,
                                connection_edge_given_level_.at(i_input_level)
                                .vec_edge_connection.at(iter_edge_index.second)
                                .vertex_connection.at(1).vertex_index));
                        }
                    } else {
                        vec_edge_others.push_back(iter_edge_index);
                    }
                }
                // keep edge connections at levels other than (i_input_level) 
                if (vec_edge_others.size() > 0) {
                    coordinate_given_level_.at(level_remove)
                        .vec_vertex_cooridinate.at(vertex_remove)
                        .vec_index_edges = std::move(vec_edge_others);
                } else {  // coordinates will be removed
                    point_remove.insert(point_index_temp);
                }
                connection_edge_given_level_.at(i_input_level)
                    .vec_edge_connection.at(edge_index_min).vertex_connection
                    = coordinate_given_level_.at(level_remove)
                    .vec_vertex_cooridinate.at(vertex_remove).parent_vertices;
                if (edge_remain.find(edge_index_min)
                    == edge_remain.end()) {
                    edge_remain.insert(edge_index_min);
                    if (dis / 2. < ds_min) {
                        edge_for_merge.insert(edge_index_min);
                    }
                }
            }

        }
    }
    // vertices whose layer levels are higher than those of the merged edges
    for (const auto& iter_point : point_remain) {
        if (point_processed.find(iter_point) == point_processed.end()) {
            point_processed.insert(iter_point);
            std::vector<std::pair<DefSizet, DefSizet>> vec_edge_others{};
            for (const auto& iter_edge_index : coordinate_given_level_
                .at(iter_point.first).vec_vertex_cooridinate
                .at(iter_point.second).vec_index_edges) {
                // only take actions for edges at (i_input_level)
                if (iter_edge_index.first == i_input_level) {
                    if (edge_remove.find(iter_edge_index.second)
                        == edge_remove.end()) {
                        edge_remove.insert(iter_edge_index.second);
                    }
                } else {
                    vec_edge_others.push_back(iter_edge_index);
                }
            }
            if (vec_edge_others.size() > 0) {
                coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.at(iter_point.second)
                    .vec_index_edges = std::move(vec_edge_others);
            } else {  // coordinates will be removed
                point_remove.insert(point_index_temp);
            }
        }
    }
    // find edges need to be removed in edge remain
    for (const auto& iter_edge : edge_remain) {
        level_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).layer_level;
        i_vertex0 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(0).vertex_index;
        level_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).layer_level;
        i_vertex1 = connection_edge_given_level_.at(i_input_level)
            .vec_edge_connection.at(iter_edge)
            .vertex_connection.at(1).vertex_index;
        if (point_processed.find(std::make_pair(level_vertex0, i_vertex0))
            != point_processed.end()
            || point_processed.find(std::make_pair(level_vertex1, i_vertex1))
            != point_processed.end()) {
            edge_remove.insert(iter_edge);
        } else {
            if (edge_for_merge.find(iter_edge) != edge_for_merge.end()) {
                ptr_surface_remain_for_merge->push_back(iter_edge);
            }
        }
    }

    // remove coordinates
    for (const auto& iter_point : point_remove) {
        i_vertex0 = coordinate_given_level_.at(iter_point.first)
            .vec_vertex_cooridinate.size();
        if (iter_point.second < i_vertex0) {
            --i_vertex0;
            while (point_remove.find(std::make_pair(iter_point.first,
                i_vertex0)) != point_remove.end() && i_vertex0 > 0) {
                // delete child relation
                level_vertex0 = coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.back().parent_vertices.at(0)
                    .layer_level;
                i_vertex0 = coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.back().parent_vertices.at(0)
                    .vertex_index;
                for (auto& iter_parent : coordinate_given_level_.at(level_vertex0)
                    .vec_vertex_cooridinate.at(i_vertex0).vec_child_vertices) {
                    if (iter_parent.first == iter_point.first
                        && iter_parent.second == iter_point.second) {
                        std::swap(iter_parent, coordinate_given_level_
                            .at(level_vertex0).vec_vertex_cooridinate
                            .at(i_vertex0).vec_child_vertices.back());
                        break;
                    }
                }
                coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.pop_back();
                --i_vertex0;
            }
            // delete coordinates and update connections
            if (iter_point.second < i_vertex0) {
                std::swap(coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.at(iter_point.second),
                    coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.at(i_vertex0));
                coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.pop_back();
                for (auto& iter_edge: coordinate_given_level_.at(iter_point.first)
                    .vec_vertex_cooridinate.at(iter_point.second).)
            }
        }
    }
}
void GeometryInfoConnection2D::DecomposeNHigerLevel(const DefSizet i_level_grid,
    const DefReal decompose_length,
    const std::unordered_map<DefSizet, bool>& map_indices_base,
    std::unordered_map<DefSizet, bool>* const ptr_map_indices_remain) {

}
//void GeometryInfoConnection2D::DecomposeOneHigerLevel(
//    const DefSizet i_level, const DefReal decompose_length) {
//    DefUint icount = 0, imax = 3;
//    // initialize coordinate and connection at current level
//    // those at the lower level
//    DefSizet vec_size = vec_coordinate_each_level_.size();
//    DefSizet layer_level = i_level - i_level_;
//    if (vec_size == layer_level) {
//        if (i_level - 1 == i_level_) {
//            vec_coordinate_each_level_.emplace_back(
//                vec_coordinate_origin_);
//            map_connection_higher_level_.emplace_back(
//                map_connection_given_level_);
//        } else {
//            vec_coordinate_each_level_.emplace_back(
//                vec_coordinate_each_level_.at(layer_level - 1));
//            map_connection_higher_level_.emplace_back(
//                map_connection_higher_level_.at(layer_level - 1));
//
//        }
//    } else if (vec_size < i_level - i_level_) {
//        io::LogWarning("Information at one level lower does not exist "
//            "for decompose geometry in function"
//            " GeometryInfoConnection::DecomposeOneHigerLevel.");
//    }
//    // initialize indices of segment will be decomposed
//    std::unordered_map<DefSizet, std::vector<DefSizet>> map_index_base;
//    DefSizet num_lines;
//    for (const auto& iter : map_connection_higher_level_.at(layer_level)) {
//        num_lines = iter.second.size();
//        std::vector<DefSizet> vec_2nd_points(num_lines);
//        for (DefSizet i = 0; i < num_lines; ++i) {
//            vec_2nd_points.at(i) = i;
//        }
//        map_index_base.insert({ iter.first, vec_2nd_points });
//    }
//
//    // splitng utill length of each segament is less than the targeted one
//    while (!map_index_base.empty()&& icount < imax) {
//        ++icount;
//        std::cout << map_index_base.size()  << " " << layer_level << "ss" << std::endl;
//        std::unordered_map<DefSizet, std::vector<DefSizet>>
//            map_index_remain;
//        DecomposeOnce(decompose_length, map_index_base,
//            &vec_coordinate_each_level_.at(layer_level),
//            &map_connection_higher_level_.at(layer_level),
//            &map_index_remain);
//        //for (const auto& iter : map_connection_higher_level_.at(layer_level)) {
//        //    for (const auto& iter_vec : iter.second) {
//        //        std::cout << iter.first << " " << iter_vec.index_second << satd
//        //    }
//      
//        //}
//        map_index_base.swap(map_index_remain);
//    }
//    if (icount >= imax) {
//        io::LogWarning("Procedure of decomposing geometry exceed"
//            "the the limit of iteration.");
//    }
//}
///**
//* @brief   function to add point between two linked vertices
//* @param[in] decompose_length  targeted length less than
//*            which a midpoint will be added
//* @param[in] map_indices_base  indices of linked vertices
//* @param[out] ptr_vec_coordi coordinates of existing and added vertices
//* @param[out] ptr_vec_connect connection of existing and added vertices
//* @param[out] ptr_map_indices_remain indices of linked vertices whose length
//*                                    is greater than the targeted one
//* @note
//*/
//void GeometryInfoConnection2D::DecomposeOnce(const DefReal decompose_length,
//    const std::unordered_map<DefSizet, std::vector<DefSizet>>&
//    map_indices_base,
//    std::vector<GeometryCoordinate>* const ptr_vec_coordi,
//    std::unordered_map<DefSizet, std::vector<GeometryConnectIndex>>*
//    const ptr_map_connect,
//    std::unordered_map<DefSizet, std::vector<DefSizet>>*
//    const ptr_map_indices_remain) {
//    io::LogWarning("The function GeometryInfoConnection::DecomposeOnce is emptty, "
//        "using derive class instead.");
//}
///**
//* @brief   function to add point between two linked vertices
//* @param[in] decompose_length  targeted length less than
//*            which a midpoint will be added
//* @param[in] map_indices_base  indices of linked vertices
//* @param[out] ptr_vec_coordi coordinates of existing and added vertices
//* @param[out] ptr_vec_connect connection of existing and added vertices
//* @param[out] ptr_map_indices_remain indices of linked vertices whose length
//*                                    is greater than the targeted one
//* @note
//*/
//void GeometryInfoConnect2DLine::DecomposeOnce(const DefReal decompose_length,
//    const std::unordered_map<DefSizet, std::vector<DefSizet>>&
//    map_indices_base,
//    std::vector<GeometryCoordinate>* const ptr_vec_coordi,
//    std::unordered_map<DefSizet, std::vector<GeometryConnectIndex>>*
//    const ptr_map_connect,
//    std::unordered_map<DefSizet, std::vector<DefSizet>>*
//    const ptr_map_indices_remain) {
//
//    // add coordinates and connections at a level higher
//    GeometryCoordinate coordinate_temp;
//    coordinate_temp.coordinate = { 0,0 };
//    GeometryConnectIndex connect_temp;
//    DefReal decompose_threshold = decompose_length * decompose_factor_;
//    DefSizet index_second, coordi_max, connect_max;
//    DefReal x_start, x_end, y_start, y_end, distance;
//    for (const auto& index_1st_point : map_indices_base) {
//        x_start = ptr_vec_coordi->at(index_1st_point.first).coordinate[kXIndex];
//        y_start = ptr_vec_coordi->at(index_1st_point.first).coordinate[kYIndex];
//        for (const auto index_2nd_point : index_1st_point.second) {
//            if (!ptr_map_connect->at(index_1st_point.first).
//                at(index_2nd_point).bool_decomposed) {
//                index_second = ptr_map_connect->
//                    at(index_1st_point.first).at(index_2nd_point).index_second;
//                x_end = ptr_vec_coordi->at(index_second).coordinate[kXIndex];
//                y_end = ptr_vec_coordi->at(index_second).coordinate[kYIndex];
//                distance = std::sqrt(Square(x_start - x_end)
//                    + Square(y_start - y_end));
//                if (distance > decompose_threshold) {
//                    // add a new point
//                    coordinate_temp.coordinate[kXIndex] = (x_start + x_end) / 2;
//                    coordinate_temp.coordinate[kYIndex] = (y_start + y_end) / 2;
//                    ptr_vec_coordi->emplace_back(coordinate_temp);
//                    // update connection information
//                    coordi_max = ptr_vec_coordi->size() - 1;
//                    ptr_map_connect->at(index_1st_point.first).
//                        at(index_2nd_point).index_second = coordi_max;
//                    connect_temp.index_second = coordi_max;
//                    if (ptr_map_connect->find(index_second) != ptr_map_connect->end()) {
//                        ptr_map_connect->at(index_second).emplace_back(connect_temp);
//                    } else {
//                        ptr_map_connect->insert({ index_second, {connect_temp} });
//                    }
//                    connect_max = ptr_map_connect->at(index_second).size() - 1;
//                    if (distance / 2. > decompose_threshold + kEps) {
//                        ptr_map_connect->at(index_1st_point.first).
//                            at(index_2nd_point).bool_decomposed = false;
//                        ptr_map_connect->at(index_second).
//                            at(connect_max).bool_decomposed = false;
//                        if (ptr_map_indices_remain->find(index_1st_point.first)
//                            == ptr_map_indices_remain->end()) {
//                            ptr_map_indices_remain->insert(
//                                { index_1st_point.first , {index_2nd_point} });
//                        } else {
//                            ptr_map_indices_remain->at(index_1st_point.first).
//                                emplace_back(index_2nd_point);
//                        }
//                        if (ptr_map_indices_remain->find(coordi_max)
//                            == ptr_map_indices_remain->end()) {
//                            ptr_map_indices_remain->insert(
//                                { index_second , {connect_max} });
//                        } else {
//                            ptr_map_indices_remain->at(index_second).
//                                emplace_back(connect_max);
//                        }
//                    }
//                    else {
//                        ptr_map_connect->at(index_1st_point.first).
//                            at(index_2nd_point).bool_decomposed = true;
//                        ptr_map_connect->at(index_second).
//                            at(connect_max).bool_decomposed = true;
//                    }
//                }
//            }
//        }
//    }
//}
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
