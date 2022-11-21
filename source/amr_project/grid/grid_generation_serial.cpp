////  Copyright (c) 2022, Zhengliang Liu
////  All rights reserved
//
///**
//* @file grid_generation_serial.cpp
//* @author Zhengliang Liu
//* @brief functions used to generate grid serially.
//* @date  2022-7-24
//* @note  functions from geometry_manager will be called.
//*/
//#include <fstream>
//#include "auxiliary_inline_func.h"
//#include "grid/grid_manager.h"
//#include "io/log_write.h"
//#ifdef ENABLE_MPI
//#include "mpi/mpi_manager.h"
//#endif  // ENABLE_MPI
//namespace rootproject {
//namespace amrproject {
///**
//* @brief   function to generate grid for all levels of refinement.
//* @param[in]  vec_geo_info vetor containing vertexers
//*             to instance storing geometry information.
//* @note
//*/
//void GridManagerInterface::GenerateGridSerial(
//    const std::vector<std::shared_ptr
//    <criterion::GeometryInfoInterface>>& vec_geo_info) {
//    // maps to store nodes at each refiement level
//    std::vector<DefMap<DefUint>> vec_map_sfbitset_at_lower_level(
//        k0MaxLevel_ + 1, DefMap<DefUint>{});
//    std::vector<DefMap<DefUint>> vec_map_tracking_node(
//        k0MaxLevel_, DefMap<DefUint>{});
//    DefMap<DefUint> map_nodes_pre_iter;
//
//    // generate tracking and ghost nodes based on geometries
//    DefSizet i_geo = 0;
//    for (auto& iter : vec_geo_info) {
//        SearchingForNodesBasedOnGeometry(k0MaxLevel_, i_geo, iter,
//            &vec_map_sfbitset_at_lower_level.at(k0MaxLevel_),
//            &vec_map_tracking_node);
//        ++i_geo;
//    }
//    // find tracking node at the highest refinement level
//    for (DefSizet i_level = k0MaxLevel_; i_level > 0; --i_level) {
//        DefSizet i_level_lower = i_level - 1;
//
//        // generate grid by searhing for nodes at the lower level. The missing
//        // nodes will be added during instantiation
//        std::vector<DefUint> vec_number_of_extend_layer_at_lower_level_neg =
//        { (k0IntExtend_.at(i_level) + k0XIntExtendNegative_.at(i_level)) / 2,
//          (k0IntExtend_.at(i_level) + k0YIntExtendNegative_.at(i_level)) / 2 };
//        std::vector<DefUint> vec_number_of_extend_layer_at_lower_level_pos =
//        { (k0IntExtend_.at(i_level) + k0XIntExtendPositive_.at(i_level)) / 2,
//          (k0IntExtend_.at(i_level) + k0YIntExtendPositive_.at(i_level)) / 2 };
//        if (k0GridDims_ == 3) {
//            vec_number_of_extend_layer_at_lower_level_neg.push_back(
//                (k0IntExtend_.at(i_level)
//                    + k0ZIntExtendNegative_.at(i_level)) / 2);
//            vec_number_of_extend_layer_at_lower_level_pos.push_back(
//                (k0IntExtend_.at(i_level)
//                    + k0ZIntExtendPositive_.at(i_level)) / 2);
//        }
//        GenerateNodeLayerByLayer(i_level - 1,
//            vec_number_of_extend_layer_at_lower_level_neg,
//            vec_number_of_extend_layer_at_lower_level_pos, map_nodes_pre_iter,
//            &vec_map_sfbitset_at_lower_level.at(i_level),
//            &vec_map_sfbitset_at_lower_level.at(i_level_lower));
//        map_nodes_pre_iter.clear();
//        map_nodes_pre_iter = vec_map_sfbitset_at_lower_level.at(i_level_lower);
//    }
//
//    for (DefSizet i_level = 1; i_level < k0MaxLevel_; ++i_level) {
//        for (auto iter = vec_map_sfbitset_at_lower_level.at(i_level).begin();
//          iter != vec_map_sfbitset_at_lower_level.at(i_level).end(); ++iter) {
//            AddNodesInstanceBasedOnLowerLevel(iter->first,
//                vec_map_sfbitset_at_lower_level.at(i_level),
//                vec_ptr_grid_info_.at(i_level));
//        }
//    }
//}
//
///**
//* @brief   function to find nodes according to the given geometry.
//* @param[in]  i_level current refinement level.
//* @param[in]  i_geo index of the geometry (only for write log).
//* @param[in]  ptr_geo_info vertexer to instance storing geometry information.
//* @param[out]  map_sfbitset_at_lower_level
//*                  nodes with bitsets at the lower level.
//* @param[out]  ptr_vec_map_tracking_nodes
//*                  the highest level bitset of tacking nodes
//                   geometry of different levels.
//* @note
//*/
//void GridManagerInterface::SearchingForNodesBasedOnGeometry(
//    const DefSizet i_level, const DefSizet i_geo,
//    const std::shared_ptr<criterion::GeometryInfoInterface> ptr_geo_info,
//    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level,
//    std::vector<DefMap<DefUint>>* const ptr_vec_map_tracking_nodes) {
//
//    // find tracking nodes
//    if (ptr_geo_info->node_type_.empty()) {
//        io::LogError("Can't find corresponding tracking node"
//            " type for geometry: "+ std::to_string(i_geo)
//            + " in SearchingForNodesBasedOnGeometry.");
//    } else if (ptr_geo_info->geometry_enclosed_type_
//        == criterion::EGeometryEnclosedType::kOpen) {
//        FindTrackingNodesBaseOnGeovertexs(i_geo,
//            i_level, ptr_geo_info->i_level_,
//            ptr_geo_info->vec_coordinate_origin_,
//            ptr_map_sfbitset_at_lower_level,
//            &(ptr_vec_map_tracking_nodes->at(ptr_geo_info->i_level_)));
//    } else {  // a temperal map is used to classify node types
//        DefMap<DefUint> map_sfbiset_for_each_gemetry;
//        FindTrackingNodesBaseOnGeovertexs(i_geo,
//            i_level, ptr_geo_info->i_level_,
//            ptr_geo_info->vec_coordinate_origin_,
//            &map_sfbiset_for_each_gemetry,
//            &(ptr_vec_map_tracking_nodes->at(ptr_geo_info->i_level_)));
//        // create layers near the geometry to classify
//        // regions of different types
//        DefMap<DefUint> map_nodes_uncolored, map_nodes_colored;
//        IdentifyTypeOfLayerByFloodFill(i_level, i_geo,
//            ptr_geo_info, map_sfbiset_for_each_gemetry,
//            &map_nodes_uncolored, &map_nodes_colored);
//        // add temperal map to map_sfbitset_at_lower_level
//        for (auto iter = map_sfbiset_for_each_gemetry.begin();
//            iter != map_sfbiset_for_each_gemetry.begin(); ++iter) {
//            if (ptr_map_sfbitset_at_lower_level->find(iter->first) ==
//                ptr_map_sfbitset_at_lower_level->end()) {
//                ptr_map_sfbitset_at_lower_level->insert(
//                    { iter->first, iter->second });
//            } else {
//                ptr_map_sfbitset_at_lower_level->at(iter->first)
//                    |= iter->second;
//            }
//        }
//    }
//}
//
///**
//* @brief   function to find tracking nodes based on geometry vertexs.
//* @param[in]  i_geo index of the geometry (only for write log).
//* @param[in]  i_level current refinement level.
//* @param[in]  i_geo_level level of geometry.
//* @param[in]  vec_geo_coordi geometry coordinates.
//* @param[out]  map_sfbitset_at_lower_level  
//*                  nodes withbit_sets at the lower level.
//* @param[out]  map_tracking_node
//*                 tracking nodes at a given level.
//*/
//void GridManagerInterface::FindTrackingNodesBaseOnGeoVertics(const DefSizet i_geo,
//    const DefSizet i_level, const DefSizet i_geo_level,
//    const std::vector<criterion::GeometryCoordinate>& vec_geo_coordi,
//    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level,
//    DefMap<DefUint>* const ptr_map_tracking_node) {
//
//    DefReal coordi_for_compare;
//    std::vector<DefReal> coordinate_min(vec_geo_coordi.begin()->coordinate);
//    std::vector<DefReal> coordinate_max(vec_geo_coordi.begin()->coordinate);
//    // set size of vector according to dimension
//    std::vector<DefSFBitset> vec_node_tracking, vec_node_tracking_lower_level;
//    if (k0GridDims_ == 2) {
//        vec_node_tracking = std::vector<DefSFBitset>(4, 0);
//    } else {
//        vec_node_tracking = std::vector<DefSFBitset>(8, 0);
//    }
//    vec_node_tracking_lower_level = vec_node_tracking;
//
//    // select dimension specified function
//    void(GridManagerInterface:: * ptrfunc_find_tracking_nodes_containing_geo_vertex)(
//        const std::vector<DefReal>&, const std::vector<DefReal>&,
//        std::vector<DefSFBitset>*) = nullptr;
//    if (k0GridDims_ == 2) {
//#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_find_tracking_nodes_containing_geo_vertex =
//            &GridManagerInterface::FindTrackingNodesContainingAGeometryCoordinate2D;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    } else if (k0GridDims_ == 3) {
//#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_find_tracking_nodes_containing_geo_vertex =
//            &GridManagerInterface::FindTrackingNodesContainingAGeometryCoordinate3D;
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//    }
//
//    // find tracking nodes based on gometry vertexs
//    DefUint node_flag = kFlag0_;
//    if (i_level > i_geo_level) {
//        node_flag = kFlagBitLowerLevel_;
//    }
//
//    for (auto iter = vec_geo_coordi.begin();
//        iter != vec_geo_coordi.end(); ++iter) {
//        // find minimum and maximum coodirnates of the geometry
//        for (DefUint i_dim = 0; i_dim < k0GridDims_; ++i_dim) {
//            coordi_for_compare = iter->coordinate.at(i_dim);
//            if (coordinate_min.at(i_dim) > coordi_for_compare) {
//                coordinate_min.at(i_dim) = coordi_for_compare;
//            } else if (coordinate_max.at(i_dim) < coordi_for_compare) {
//                coordinate_max.at(i_dim) = coordi_for_compare;
//            }
//        }
//
//        // find tracking nodes when grid and geometry are at the same level
//        if (i_level == i_geo_level) {
//            (this->*ptrfunc_find_tracking_nodes_containing_geo_vertex)(
//                vec_ptr_grid_info_.at(i_level)->grid_space_, (*iter).coordinate,
//                &vec_node_tracking);
//            for (DefSFBitset& iter_sfbitset : vec_node_tracking) {
//                if (ptr_map_tracking_node->find(iter_sfbitset)
//                    == ptr_map_tracking_node->end()) {
//                    ptr_map_tracking_node->insert(
//                        { iter_sfbitset, 1 });
//                }
//            }
//        }
//
//        // find tracking nodes at lower level
//        (this->*ptrfunc_find_tracking_nodes_containing_geo_vertex)(
//            vec_ptr_grid_info_.at(i_level - 1)->grid_space_,
//            (*iter).coordinate, &vec_node_tracking_lower_level);
//        for (const auto& iter_sfbitset : vec_node_tracking_lower_level) {
//            if (ptr_map_sfbitset_at_lower_level->find(iter_sfbitset)
//                == ptr_map_sfbitset_at_lower_level->end()) {
//                ptr_map_sfbitset_at_lower_level->insert(
//                    { iter_sfbitset, node_flag });
//            } else {
//                ptr_map_sfbitset_at_lower_level->at(iter_sfbitset)
//                    |= node_flag;
//            }
//        }
//    }
//
//    int status = CheckIfvertexOutsideDomain(coordinate_min, coordinate_max,
//        k0RealOffest_, k0DomainSize_);
//
//    if (status == 1) {
//        io::LogError("The minimum coordinates of Geometry "
//            + std::to_string(i_geo) + " is less than the domain");
//    } else if (status == 2) {
//        io::LogError("The maximum coordinates of Geometry "
//            + std::to_string(i_geo) + " is less than the domain");
//    }
//}
//
///**
//* @brief   function to add layers near geometry vertexs.
//* @param[in] i_level refinement level
//* @param[in]  i_geo index of the geometry (only for write log).
//* @param[in]  ptr_geo_info vertexer to instance storing geometry information.
//* @param[in]  map_sfbitset   existing nodes.
//* @param[out] ptr_map_nodes_uncolored   nodes haven't been colored.
//* @param[out] ptr_map_nodes_colored   nodes havebeen colored.
//*/
//void GridManagerInterface::IdentifyTypeOfLayerByFloodFill(
//    const DefSizet i_level, const DefSizet i_geo,
//    const std::shared_ptr<criterion::GeometryInfoInterface> ptr_geo_info,
//    const DefMap<DefUint>&  map_sfbitset,
//    DefMap<DefUint>* const ptr_map_nodes_uncolored,
//    DefMap<DefUint>* const ptr_map_nodes_colored) {
//
//    DefSizet i_geo_level = ptr_geo_info->i_level_;
//
//    // identify nodes layer type through flood fill,
//    // which is conducted on ptr_map_sfbitset_at_lower_level
//    // where the refinement level is i_level - 1.
//    std::vector<DefReal> flood_fill_origin;
//    // set cooridate for finding flood fill starting vertex the same
//    // as geometry centor if flood_fill_origin_ is not given.
//    if (flood_fill_origin.empty()) {
//        flood_fill_origin = ptr_geo_info->geometry_center_;
//    } else {
//        flood_fill_origin = ptr_geo_info->flood_fill_origin_;
//    }
//
//    // idengtify nodes type using flood fill method.
//    DefSFBitset sfbitset_start;
//    if (k0GridDims_ == 2) {
//        sfbitset_start = FindStartingvertexForFloodFill2D(
//            i_level, i_geo, flood_fill_origin, map_sfbitset);
//    } else if (k0GridDims_ == 3) {
//        sfbitset_start = FindStartingvertexForFloodFill3D(
//            i_level, i_geo, flood_fill_origin, map_sfbitset);
//    }
//    FloodFillOneLayer(sfbitset_start, map_sfbitset,
//         ptr_map_nodes_uncolored, ptr_map_nodes_colored);
//}
///**
//* @brief   function to add k0IntGhostExtend_ layers.
//* @param[in]  i_level   level of refinement.
//* @param[in]  map_tracking_node_temp
//*                  tracking nodes classified to find ghost nodes.
//* @param[out]  map_sfbitset_at_lower_level
//*                  nodes withbit_sets at the lower level.
//* @note  add the number of layer around the geometry in all the directions.
//*        Used to set number of layers inside the geometry different from that
//*        outsite the geometry.
//*/
//void GridManagerInterface::ExtendSameNumberOfLayer(const DefSizet i_level,
//    const DefMap<DefUint>& map_tracking_node_temp,
//    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level) {
//
//    DefSizet i_level_lower = i_level - 1;
//    DefSFBitset sfbitset_at_lower_level;
//    // set extended layer the same in all directions
//    std::vector<DefLUint> vec_extend_neg
//    (k0GridDims_, k0IntInnerExtend_.at(i_level));
//    std::vector<DefLUint> vec_extend_pos
//    (k0GridDims_, k0IntInnerExtend_.at(i_level));
//
//    if (k0GridDims_ == 2) {
//#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
//        DefSFBitset sfbitset_temp_x, sfbitset_temp_y;
//        for (auto iter = map_tracking_node_temp.begin();
//            iter != map_tracking_node_temp.end(); ++iter) {
//            // sfbitset at the lower finement level
//            sfbitset_at_lower_level =
//                SFBitsetToNLowerLevel2D(1, iter->first);
//            ResetExtendLayerBasedOnDomainSize2D(
//                i_level_lower, sfbitset_at_lower_level,
//                &vec_extend_neg, &vec_extend_pos);
//
//            // find sfbitset in a region
//            std::array<DefLUint,2> indices(k0GridDims_, 0);
//            SFBitsetComputeIndices2D(sfbitset_at_lower_level, &indices);
//            std::array<DefLUint, 2> indices_reset;
//            indices_reset = { indices[kXIndex] - vec_extend_neg[kXIndex],
//            indices[kYIndex] - vec_extend_neg[kYIndex] };
//            DefSFBitset sfbitset_y = SFBitsetEncoding2D(indices_reset);
//            DefSFBitset sfbitset_x;
//            for (DefLUint iy = 0; iy <= vec_extend_neg[kYIndex]
//                + vec_extend_pos[kYIndex]; ++iy) {
//                sfbitset_x = sfbitset_y;
//                for (DefLUint ix = 0; ix <= vec_extend_neg[kXIndex]
//                    + vec_extend_pos[kXIndex]; ++ix) {
//                    if (ptr_map_sfbitset_at_lower_level->find(sfbitset_x)
//                        == ptr_map_sfbitset_at_lower_level->end()) {
//                        ptr_map_sfbitset_at_lower_level->insert(
//                            { sfbitset_x, kFlag0_ });
//                    }
//                    sfbitset_x = FindXPos2D(sfbitset_x);
//                }
//                sfbitset_y = FindYPos2D(sfbitset_y);
//            }
//        }
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    }  else if (k0GridDims_ == 3) {
//#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
//        DefSFBitset sfbitset_temp_x, sfbitset_temp_y, sfbitset_temp_z;
//        for (auto iter = map_tracking_node_temp.begin();
//            iter != map_tracking_node_temp.end(); ++iter) {
//            // sfbitset at the lower finement level
//            sfbitset_at_lower_level =
//                SFBitsetToNLowerLevel3D(1, iter->first);
//            ResetExtendLayerBasedOnDomainSize3D(
//                i_level_lower, sfbitset_at_lower_level,
//                &vec_extend_neg, &vec_extend_pos);
//
//            // find sfbitset in a region
//            std::array<DefLUint, 3> indices(k0GridDims_, 0);
//            SFBitsetComputeIndices3D(sfbitset_at_lower_level, &indices);
//            DefSFBitset sfbitset_z = SFBitsetEncoding3D(
//                { indices[kXIndex] - vec_extend_neg[kXIndex],
//                indices[kYIndex] - vec_extend_neg[kYIndex],
//                indices[kZIndex] - vec_extend_neg[kZIndex] });
//            DefSFBitset sfbitset_x, sfbitset_y;
//            for (DefLUint iz = 0; iz <= vec_extend_neg[kZIndex]
//                + vec_extend_pos[kZIndex]; ++iz) {
//                sfbitset_y = sfbitset_z;
//                for (DefLUint iy = 0; iy <= vec_extend_neg[kYIndex]
//                    + vec_extend_pos[kYIndex]; ++iy) {
//                    sfbitset_x = sfbitset_y;
//                    for (DefLUint ix = 0; ix <= vec_extend_neg[kXIndex]
//                        + vec_extend_pos[kXIndex]; ++ix) {
//                        if (ptr_map_sfbitset_at_lower_level->find(sfbitset_x)
//                            == ptr_map_sfbitset_at_lower_level->end()) {
//                            ptr_map_sfbitset_at_lower_level->insert(
//                                { sfbitset_x, kFlag0_ });
//                        }
//                        sfbitset_x = FindXPos3D(sfbitset_x);
//                    }
//                    sfbitset_y = FindYPos3D(sfbitset_y);
//                }
//                sfbitset_z = FindZPos3D(sfbitset_z);
//            }
//        }
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//    }
//}
///**
//* @brief   function to do flood fill on existing nodes around geometry vertexs
//*          (only two layers)
//* @param[in]  sfbitset_start   sfbitset corresponding to the starting vertex.
//* @param[in]  map_nodes_exist nodes exist for flood fill.
//* @param[out]  ptr_map_nodes_uncolored
//*                  nodes not colored by the flood fill method.
//* @param[out]  ptr_map_nodes_colored
//*                  nodes colored by the flood fill method.
//* @node ptr_map_nodes_colored contains nodes on the inner layer, while
//*       ptr_map_nodes_uncolored contains nodes on the outer layer.
//*       The diagonal nodes are considered on the outer layer0 
//*/
////  x   x           //  o is geometry vertexs
////    o             //  x are nodes of containing the geometry vertexs
////  x  (x)  x       //  (x) is the diagonal node
////        o
////      x   x
//void GridManagerInterface::FloodFillOneLayer(const DefSFBitset& sfbitset_start,
//    const DefMap<DefUint>&  map_nodes_exist,
//    DefMap<DefUint>* const ptr_map_nodes_uncolored,
//    DefMap<DefUint>* const ptr_map_nodes_colored) {
//
//    // select dimesion specified function
//    void(GridManagerInterface:: * ptrfunc_push_back_sfbitset)(
//        const DefSFBitset & sfbitset_in,
//        std::vector<DefSFBitset>*ptr_vec_code_stk) = nullptr;
//    if (k0GridDims_ == 2) {
//#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_push_back_sfbitset =
//            &GridManagerInterface::PushBackSFBitsetInFloodFill2D;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//        } else if (k0GridDims_ == 3) {
//#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_push_back_sfbitset =
//            &GridManagerInterface::PushBackSFBitsetInFloodFill3D;
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//    }
//
//    // create a copy of map_nodes_exist and add one layer for flood fill
//    DefUint flag_bit_ghost = 1;
//    DefUint flag_bit_colored = 1 << 3;
//    DefUint flag_bit_exist = flag_bit_colored | (1 << 4);
//    DefMap<DefUint> map_nodes_temp(map_nodes_exist);
//    std::vector<DefSFBitset> vec_bitset;
//
//    // choose dimesion specified functions
//    void (*bitset_find_all_neighbours_ptr)(const DefSFBitset&,
//        std::vector<DefSFBitset>*) = nullptr;
//    if (k0GridDims_ == 2) {
//        vec_bitset = std::vector<DefSFBitset>(9, 0);
//    } else {
//        vec_bitset = std::vector<DefSFBitset>(27, 0);
//    }
//    // generate one more layer around the tracking nodes
//    for (auto iter = map_nodes_exist.begin();
//        iter != map_nodes_exist.end(); ++iter) {
//        map_nodes_temp.at(iter->first) = flag_bit_exist;
//        (*bitset_find_all_neighbours_ptr)(iter->first, &vec_bitset);
//        for (const auto& iter_neigbour : vec_bitset) {
//            if (map_nodes_temp.find(iter_neigbour)
//                == map_nodes_temp.end()) {
//                map_nodes_temp.insert({ iter_neigbour, flag_bit_ghost });
//            }
//        }
//    }
//
//    // flood fill
//    std::vector<DefSFBitset> vec_sfbitset_stk;
//    DefUint i = 0;
//    DefSFBitset sfbitset_seed;
//
//    vec_sfbitset_stk.push_back(sfbitset_start);
//    while (!vec_sfbitset_stk.empty() && i < imax_flood_fill) {
//        sfbitset_seed = vec_sfbitset_stk.back();
//        vec_sfbitset_stk.pop_back();
//        ++i;
//        if (map_nodes_temp.find(sfbitset_seed) == map_nodes_temp.end()) {
//        } else if ((map_nodes_temp.at(sfbitset_seed) & flag_bit_exist)
//            == flag_bit_exist) {
//            ptr_map_nodes_colored->insert({ sfbitset_seed, kFlag0_ });
//        } else if ((map_nodes_temp.at(sfbitset_seed) & flag_bit_colored)
//            != flag_bit_colored) {
//            // color the node
//            map_nodes_temp.at(sfbitset_seed) |= flag_bit_colored;
//            // add neighbouring nodes to seed
//            (this->*ptrfunc_push_back_sfbitset)(
//                sfbitset_seed, &vec_sfbitset_stk);
//        }
//    }
//
//    if (i == imax_flood_fill) {
//        io::LogInfo("Iteration of flood fill exceed preset limits.");
//    }
//
//    // map_nodes_uncolored = map_nodes_exist - map_nodes_colored
//    for (auto iter = map_nodes_exist.begin();
//        iter != map_nodes_exist.end(); ++iter) {
//        if (ptr_map_nodes_colored->find(iter->first)
//            == ptr_map_nodes_colored->end()) {
//            ptr_map_nodes_uncolored->insert({ iter->first, kFlag0_ });
//        }
//    }
//}
///**
//* @brief   function to generate nodes extened from the given layer
//* @param[in]  i_level   level of refinement.
//* @param[in]  vec_number_of_extend_layer_neg   number of layer need 
//                           to be extened in negative direction.
//* @param[in]  vec_number_of_extend_layer_pos   number of layer need
//                           to be extened in positive direction.
//* @param[out]  ptr_map_nodes_uncolored
//*                  nodes not colored by the flood fill method.
//* @param[out]  ptr_map_nodes_colored
//*                  nodes colored by the flood fill method.
//*/
//void GridManagerInterface::GenerateNodeLayerByLayer(const DefSizet i_level,
//    const std::vector<DefUint>& vec_number_of_extend_layer_neg,
//    const std::vector<DefUint>& vec_number_of_extend_layer_pos,
//    const DefMap<DefUint>& map_nodes_starting_layer,
//    DefMap<DefUint>* const ptr_map_nodes_exist,
//    DefMap<DefUint>* const map_two_outmost_layer_at_lower_level) {
//
////    // identify in which direction the grid is expanded the most
////    DefUint extend_max = vec_number_of_extend_layer_neg[kXIndex];
////    DefUint extend_min = vec_number_of_extend_layer_neg[kXIndex];
////    if (vec_number_of_extend_layer_pos[kXIndex] > extend_max) {
////        extend_max = vec_number_of_extend_layer_pos[kXIndex];
////    } else if (vec_number_of_extend_layer_pos[kXIndex] < extend_min) {
////        extend_min = vec_number_of_extend_layer_pos[kXIndex];
////    }
////    if (vec_number_of_extend_layer_neg[kYIndex] > extend_max) {
////        extend_max = vec_number_of_extend_layer_neg[kYIndex];
////    } else if (vec_number_of_extend_layer_neg[kYIndex] < extend_min) {
////        extend_min = vec_number_of_extend_layer_neg[kYIndex];
////    }
////    if (vec_number_of_extend_layer_pos[kYIndex] > extend_max) {
////        extend_max = vec_number_of_extend_layer_pos[kYIndex];
////    } else if (vec_number_of_extend_layer_pos[kYIndex] < extend_min) {
////        extend_min = vec_number_of_extend_layer_pos[kYIndex];
////    }
////    if (k0GridDims_ == 3) {
////        if (vec_number_of_extend_layer_neg[kZIndex] > extend_max) {
////            extend_max = vec_number_of_extend_layer_neg[kZIndex];
////        } else if (vec_number_of_extend_layer_neg[kZIndex] < extend_min) {
////            extend_min = vec_number_of_extend_layer_neg[kZIndex];
////        }
////        if (vec_number_of_extend_layer_pos[kZIndex] > extend_max) {
////            extend_max = vec_number_of_extend_layer_pos[kZIndex];
////        } else if (vec_number_of_extend_layer_pos[kZIndex] < extend_min) {
////            extend_min = vec_number_of_extend_layer_pos[kZIndex];
////        }
////    }
////
////    // select dimension specified functions
////    void(SFBitsetAux:: *ptrfunc_check_if_node_not_on_domain_boundary)(
////        const DefSFBitset & sfbitset_in,
////        const std::vector<DefSFBitset>&vec_sfbitset_min,
////        const std::vector<DefSFBitset>&vec_sfbitset_max,
////        std::vector<bool>*ptr_vec_bool_not_at_boundary_neg,
////        std::vector<bool>*ptr_vec_bool_not_at_boundary_pos) = nullptr;
////    bool(GridManagerInterface:: * ptrfunc_add_nodes_around_a_given_node)(
////        const DefSFBitset & sfbitset_in,
////        const std::vector<bool>&vec_bool_need_add_neg,
////        const std::vector<bool>&vec_bool_need_add_pos,
////        DefMap<DefLUint>*ptr_map_nodes_exist,
////        DefMap<DefLUint>*ptr_map_nodes_aded) = nullptr;
////    std::vector<DefSFBitset> vec_sfbitset_min, vec_sfbitset_max;
////    SFBitsetToNHigherLevel(i_level,
////        sfbitset_aux_.k0SFBitsetMin,  &vec_sfbitset_min);
////    SFBitsetToNHigherLevel(i_level,
////        sfbitset_aux_.k0SFBitsetMin, &vec_sfbitset_min);
////
////    if (k0GridDims_ == 2) {
////#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
////        ptrfunc_check_if_node_not_on_domain_boundary =
////            &SFBitsetAux::SFBitsetNotOnDomainBoundary2D;
////        ptrfunc_add_nodes_around_a_given_node =
////            &GridManagerInterface::AddNodesAroundAGivenNode2D;
////#endif  // DEBUG_DISABLE_2D_FUNCTIONS
////    } else if (k0GridDims_ == 3) {
////#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
////        ptrfunc_check_if_node_not_on_domain_boundary =
////            &SFBitsetAux::SFBitsetNotOnDomainBoundary3D;
////        ptrfunc_add_nodes_around_a_given_node =
////            &GridManagerInterface::AddNodesAroundAGivenNode3D;
////#endif  // DEBUG_DISABLE_3D_FUNCTIONS
////    }
////
////    // extend mesh layer based on previous layer
////    std::vector<bool> vec_bool_extend_neg(k0GridDims_, true),
////        vec_bool_extend_pos(k0GridDims_, true);
////    DefMap<DefUint> map_outmost_layer(map_nodes_starting_layer);
////    DefMap<DefUint> map_nodes_on_any_boundaries{};
////    for (unsigned int i_extend = 1; i_extend < extend_max + 1; ++i_extend) {
////        if (i_extend > extend_min) {
////            if (i_extend > vec_number_of_extend_layer_neg[kXIndex]) {
////                vec_bool_extend_neg.at(kXIndex) = false;
////            }
////            if (i_extend > vec_number_of_extend_layer_pos[kXIndex]) {
////                vec_bool_extend_pos.at(kXIndex) = false;
////            }
////            if (i_extend > vec_number_of_extend_layer_neg[kYIndex]) {
////                vec_bool_extend_neg.at(kYIndex) = false;
////            }
////            if (i_extend > vec_number_of_extend_layer_pos[kYIndex]) {
////                vec_bool_extend_pos.at(kYIndex) = false;
////            }
////            if (k0GridDims_ == 3) {
////                if (i_extend > vec_number_of_extend_layer_neg[kZIndex]) {
////                    vec_bool_extend_neg.at(kZIndex) = false;
////                }
////                if (i_extend > vec_number_of_extend_layer_pos[kZIndex]) {
////                    vec_bool_extend_pos.at(kZIndex) = false;
////                }
////            }
////        }
////        DefMap<DefUint> map_current_layer(map_outmost_layer);
////        map_outmost_layer.clear();
////        for (auto iter = map_current_layer.begin();
////            iter != map_current_layer.end(); ++iter) {
////            // check if node is not on the domain boundary
////            (sfbitset_aux_.*ptrfunc_check_if_node_not_on_domain_boundary)(
////                iter->first, vec_sfbitset_min, vec_sfbitset_max,
////                &vec_bool_extend_neg, &vec_bool_extend_pos);
////            // add nodes in directions specified by vec_bool_extend_neg
////            // and vec_bool_extend_pos. Return true if there exist directions
////            // in which nodes do not need to be added
////            if ((this->*ptrfunc_add_nodes_around_a_given_node)(iter->first,
////                vec_bool_extend_neg, vec_bool_extend_pos,
////                ptr_map_nodes_exist, &map_outmost_layer)) {
////                map_nodes_on_any_boundaries.insert({ iter->first, kFlag0_ });
////            }
////        }
////    }
//}
//#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
///**
//* @brief   function to find starting vertex for flood fill (2D)
//* @param[in]  i_level   level of refinement.
//* @param[in]  i_geo     index of the geometry (log write only).
//* @param[in]  vec_origin origin coorfinates to search node for flood fill.
//* @param[in]  map_nodes_exist nodes exist for flood fill.
//* @return     sfbitset corresponding to the starting vertex.
//*/
//DefSFBitset GridManagerInterface::FindStartingvertexForFloodFill2D(
//    const DefSizet i_level, const DefSizet i_geo,
//    const std::vector<DefReal>& vec_origin,
//    const DefMap<DefUint>& map_nodes_exist) {
//#ifdef DEBUG_CHECK_GRID
//    if (vec_origin.size() != 2) {
//        io::LogError("The dimension of vec_origin should be 2 rather than "
//            + std::to_string(vec_origin.size()) + " for geometry: "
//            + std::to_string(i_geo) + " in FindStartingvertexForFloodFill2D.");
//    }
//#endif  // DEBUG_CHECK_GRID
//    // calculate bounds of searching step based on domain boundary
//    DefLUint scale_i_level = TwoPowerN(static_cast<DefLUint>(i_level));
//    DefLUint x_index = static_cast<DefLUint>(vec_origin.at(kXIndex)
//        / (k0DomainDx_[kXIndex] / scale_i_level) + kEps);
//    DefLUint y_index = static_cast<DefLUint>(vec_origin.at(kYIndex)
//        / (k0DomainDx_[kYIndex] / scale_i_level) + kEps);
//    DefUint x_index_max =
//        k0MaxIndexOfBackgroundNode_[kXIndex] * scale_i_level - x_index;
//    DefUint y_index_max =
//        k0MaxIndexOfBackgroundNode_[kYIndex] * scale_i_level - y_index;
//
//    bool bool_find_node_for_flood_fill = false;
//    DefUint i_count = 0, count_sum = 0;;
//    DefSFBitset sfbitset_origin_vertex =
//        SFBitsetEncoding2D(std::array<DefLUint, 2>({ x_index , y_index }));
//    // search in -x direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    DefSFBitset sfbitset_temp = sfbitset_origin_vertex, sfbitset_start_vertex;
//    while (i_count < x_index) {
//       sfbitset_temp = FindXNeg2D(sfbitset_temp);
//        if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//            sfbitset_start_vertex = FindXPos2D(sfbitset_temp);
//            bool_find_node_for_flood_fill = true;
//            break;
//        }
//        ++i_count;
//    }
//    count_sum += i_count;
//    // search in -y direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//        sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < y_index) {
//           sfbitset_temp = FindYNeg2D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindYPos2D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in +x direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < x_index_max) {
//           sfbitset_temp = FindXPos2D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindXNeg2D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in +y direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < y_index_max) {
//           sfbitset_temp = FindYPos2D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindYNeg2D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//
//    if (bool_find_node_for_flood_fill) {
//        return sfbitset_start_vertex;
//    } else {
//        io::LogError("Can't find starting node for food fill in geometry: "
//            + std::to_string(i_geo) + " in FindStartingvertexForFloodFill2D"
//            + " after " + std::to_string(count_sum) + " iterations.");
//        return 1;
//    }
//}
///**
//* @brief   function to pushbit_sets into stack (2D)
//* @param[in] sfbitset_in   bitset of spacing filling code
//*                          corresponding to node has been colored.
//* @param[out]  ptr_vec_stk     stack to store seeds.
//*/
//void GridManagerInterface::PushBackSFBitsetInFloodFill2D(const DefSFBitset& sfbitset_in,
//    std::vector<DefSFBitset>* const ptr_vec_stk) {
//    // add neighbouring nodes to seed
//    ptr_vec_stk->push_back(FindXNeg2D(sfbitset_in));
//    ptr_vec_stk->push_back(FindXPos2D(sfbitset_in));
//    ptr_vec_stk->push_back(FindYNeg2D(sfbitset_in));
//    ptr_vec_stk->push_back(FindYPos2D(sfbitset_in));
//}
///**
//* @brief function to add nodes abourd a node with specified directions (2D)
//* @param[in]  sfbitset_in  bitset of a given node.
//* @param[in]  vec_bool_need_add_neg   identifier indicate
//*                whether nodes need to be added in negative directions.
//* @param[in]  vec_bool_need_adde_pos   identifier indicate
//*                whether nodes need to be added in positive directions.
//* @param[out] ptr_map_nodes_exist nodes already exist or added in this function
//* @param[out] ptr_map_nodes_added nodes are added in this function
//* @return     if exists directions in which node does not need to be added 
//*/
//bool GridManagerInterface::AddNodesAroundAGivenNode2D(
//    const DefSFBitset& sfbitset_in,
//    const std::vector<bool>& vec_bool_need_add_neg,
//    const std::vector<bool>& vec_bool_need_add_pos,
//    DefMap<DefLUint>* const ptr_map_nodes_exist,
//    DefMap<DefLUint>* const ptr_map_nodes_added) {
//    bool bool_node_not_added_in_any_diretion = false;
//    DefSFBitset sfbitset_temp0, sfbitset_temp1, sfbitset_temp2;
//    // node at (-x, 0)
//    if (vec_bool_need_add_neg[kXIndex]) {
//        sfbitset_temp0 = FindXNeg3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (-x, -y)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (-x, +y)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (+x, 0)
//    if (vec_bool_need_add_pos[kXIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (+x, -y)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (+x, +y)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, -y)
//    if (vec_bool_need_add_neg[kYIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, +y)
//    if (vec_bool_need_add_pos[kYIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    return bool_node_not_added_in_any_diretion;
//}
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
///**
//* @brief   function to do flood fill on existing nodes (3D)
//* @param[in]  i_level   level of refinement.
//* @param[in]  i_geo     index of the geometry.
//* @param[in]  vec_origin origin coorfinates to search node for flood fill.
//* @param[in]  map_nodes_exist nodes exist for flood fill.
//* @return     bit_set corresponding to the starting vertex.
//*/
//DefSFBitset GridManagerInterface::FindStartingvertexForFloodFill3D(
//    const DefSizet i_level, const DefSizet i_geo,
//    const std::vector<DefReal>& vec_origin,
//    const DefMap<DefUint>& map_nodes_exist) {
//#ifdef DEBUG_CHECK_GRID
//    if (vec_origin.size() != 3) {
//        io::LogError("The dimension of vec_origin should be 3 rather than "
//            + std::to_string(vec_origin.size()) + "for geometry: "
//            + std::to_string(i_geo) + " in FindStartingvertexForFloodFill3D.");
//    }
//#endif  // DEBUG_CHECK_GRID
//    // calculate bounds of searching step based on domain boundary
//    DefLUint scale_i_level = TwoPowerN(static_cast<DefLUint>(i_level));
//    DefLUint x_index = static_cast<DefLUint>(vec_origin.at(kXIndex)
//        / (k0DomainDx_[kXIndex] / scale_i_level) + kEps);
//    DefLUint y_index = static_cast<DefLUint>(vec_origin.at(kYIndex)
//        / (k0DomainDx_[kYIndex] / scale_i_level) + kEps);
//    DefLUint z_index = static_cast<DefLUint>(vec_origin.at(kZIndex)
//        / (k0DomainDx_[kZIndex] / scale_i_level) + kEps);
//    DefUint x_index_max =
//        k0MaxIndexOfBackgroundNode_[kXIndex] * scale_i_level - x_index;
//    DefUint y_index_max =
//        k0MaxIndexOfBackgroundNode_[kYIndex] * scale_i_level - y_index;
//    DefUint z_index_max =
//        k0MaxIndexOfBackgroundNode_[kZIndex] * scale_i_level - z_index;
//
//    bool bool_find_node_for_flood_fill = false;
//    DefUint i_count = 0, count_sum = 0;;
//    DefSFBitset sfbitset_origin_vertex = SFBitsetEncoding3D(
//        { x_index , y_index, z_index });
//    // search in -x direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    DefSFBitset sfbitset_temp = sfbitset_origin_vertex, sfbitset_start_vertex;
//    while (i_count < x_index) {
//       sfbitset_temp = FindXNeg3D(sfbitset_temp);
//       if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//           sfbitset_start_vertex = FindXPos3D(sfbitset_temp);
//           bool_find_node_for_flood_fill = true;
//           break;
//       }
//        ++i_count;
//    }
//    count_sum += i_count;
//    // search in -y direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < y_index) {
//           sfbitset_temp = FindYNeg3D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindYPos3D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in -z direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < z_index) {
//           sfbitset_temp = FindZNeg3D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindZPos3D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in +x direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < x_index_max) {
//           sfbitset_temp = FindXPos3D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindXNeg3D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in +y direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < y_index_max) {
//           sfbitset_temp = FindYPos3D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindYNeg3D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//    // search in +z direction from the vec_origin untill meet the first vertex
//    // in map_nodes_exist
//    if (!bool_find_node_for_flood_fill) {
//       sfbitset_temp = sfbitset_origin_vertex;
//        i_count = 0;
//        while (i_count < z_index_max) {
//           sfbitset_temp = FindZPos3D(sfbitset_temp);
//            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
//               sfbitset_start_vertex = FindZNeg3D(sfbitset_temp);
//                bool_find_node_for_flood_fill = true;
//                break;
//            }
//            ++i_count;
//        }
//    }
//    count_sum += i_count;
//
//    if (bool_find_node_for_flood_fill) {
//        return sfbitset_start_vertex;
//    } else {
//        io::LogError("Can't find starting node for food fill in geometry: "
//            + std::to_string(i_geo) + " in FindStartingvertexForFloodFill3D"
//            + " after " + std::to_string(count_sum) + " iterations.");
//        return 1;
//    }
//}
///**
//* @brief   function to pushbit_sets into stack (3D)
//* @param[in]   sfbitset_in   bitset of spacing filling code
//*                          corresponding to node has been colored.
//* @param[out]  ptr_vec_stk     stack to store seeds.
//*/
//void GridManagerInterface::PushBackSFBitsetInFloodFill3D(const DefSFBitset& sfbitset_in,
//    std::vector<DefSFBitset>* const ptr_vec_stk) {
//    // add neighbouring nodes to seed
//    ptr_vec_stk->push_back(FindXNeg3D(sfbitset_in));
//    ptr_vec_stk->push_back(FindXPos3D(sfbitset_in));
//    ptr_vec_stk->push_back(FindYNeg3D(sfbitset_in));
//    ptr_vec_stk->push_back(FindYPos3D(sfbitset_in));
//    ptr_vec_stk->push_back(FindZNeg3D(sfbitset_in));
//    ptr_vec_stk->push_back(FindZPos3D(sfbitset_in));
//}
///**
//* @brief function to add nodes abourd a node with specified directions (3D)
//* @param[in]  sfbitset_in  bitset of a given node.
//* @param[in]  vec_bool_need_add_neg   identifier indicate 
//*                whether nodes need to be added in negative directions.
//* @param[in]  vec_bool_need_adde_pos   identifier indicate 
//*                whether nodes need to be added in positive directions.
//* @param[out] ptr_map_nodes_exist nodes already exist or added in this function
//* @param[out] ptr_map_nodes_added nodes are added in this function
//* @return     if exists directions in which node does not need to be added 
//*/
//bool GridManagerInterface::AddNodesAroundAGivenNode3D(
//    const DefSFBitset& sfbitset_in,
//    const std::vector<bool>& vec_bool_need_add_neg,
//    const std::vector<bool>& vec_bool_need_add_pos,
//    DefMap<DefLUint>* const ptr_map_nodes_exist,
//    DefMap<DefLUint>* const ptr_map_nodes_added) {
//    bool bool_node_not_added_in_any_diretion = false;
//    DefSFBitset sfbitset_temp0, sfbitset_temp1, sfbitset_temp2;
//    // node at (-x, 0, 0)
//    if (vec_bool_need_add_neg[kXIndex]) {
//        sfbitset_temp0 = FindXNeg3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (-x, -y, 0)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (-x, +y, 0)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (+x, 0, 0)
//    if (vec_bool_need_add_pos[kXIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (+x, -y, 0)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (+x, +y, 0)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, -y, 0)
//    if (vec_bool_need_add_neg[kYIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, +y, 0)
//    if (vec_bool_need_add_pos[kYIndex]) {
//        sfbitset_temp0 = FindXPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, 0, -z)
//    if (vec_bool_need_add_neg[kZIndex]) {
//        sfbitset_temp0 = FindZNeg3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (-x, 0, -z)
//        if (vec_bool_need_add_neg[kXIndex]) {
//            sfbitset_temp1 = FindXNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (-x, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp2 = FindYNeg3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (-x, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp2 = FindYPos3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (+x, 0, -z)
//        if (vec_bool_need_add_pos[kXIndex]) {
//            sfbitset_temp1 = FindXPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (+x, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp2 = FindYNeg3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (+x, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp2 = FindYPos3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (0, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        //  node at (0, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    // node at (0, 0, +z)
//    if (vec_bool_need_add_pos[kZIndex]) {
//        sfbitset_temp0 = FindZPos3D(sfbitset_in);
//        if (ptr_map_nodes_exist->find(sfbitset_temp0)
//            == ptr_map_nodes_exist->end()) {
//            ptr_map_nodes_exist->insert({ sfbitset_temp0, kFlag0_ });
//            ptr_map_nodes_added->insert({ sfbitset_temp0, kFlag0_ });
//        }
//        // node at (-x, 0, -z)
//        if (vec_bool_need_add_neg[kXIndex]) {
//            sfbitset_temp1 = FindXNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (-x, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp2 = FindYNeg3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (-x, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp2 = FindYPos3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (+x, 0, -z)
//        if (vec_bool_need_add_pos[kXIndex]) {
//            sfbitset_temp1 = FindXPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        // node at (+x, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp2 = FindYNeg3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (+x, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp2 = FindYPos3D(sfbitset_temp1);
//            if (ptr_map_nodes_exist->find(sfbitset_temp2)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp2, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp2, kFlag0_ });
//            }
//        }
//        // node at (0, -y, -z)
//        if (vec_bool_need_add_neg[kYIndex]) {
//            sfbitset_temp1 = FindYNeg3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//        //  node at (0, +y, -z)
//        if (vec_bool_need_add_pos[kYIndex]) {
//            sfbitset_temp1 = FindYPos3D(sfbitset_temp0);
//            if (ptr_map_nodes_exist->find(sfbitset_temp1)
//                == ptr_map_nodes_exist->end()) {
//                ptr_map_nodes_exist->insert({ sfbitset_temp1, kFlag0_ });
//                ptr_map_nodes_added->insert({ sfbitset_temp1, kFlag0_ });
//            }
//        }
//    } else {
//        bool_node_not_added_in_any_diretion = true;
//    }
//    return bool_node_not_added_in_any_diretion;
//}
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//}  // end namespace amrproject
//}  // end namespace rootproject
