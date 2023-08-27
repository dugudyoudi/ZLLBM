//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_manager2d.cpp
* @author Zhengliang Liu
* @date  2022-5-23
* @brief functions used to manage grid related processes when mpi is enabled.
*/
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#ifdef ENABLE_MPI
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
namespace rootproject {
namespace amrproject {
/**
 * @brief function to get spacing fill codes whose bits are 1 for the given level and 0 for the background
 * @param[in] i_level given refinement level
 * @param[in] bitset_aux3d class manage space filling curves.
 * @param[out] ptr_last_ones pointer to spacing fill codes
 * @throws ErrorType if the size of last_ones is not 3
 */
void MpiManager::GetNLevelCorrespondingOnes3D(const DefAmrIndexUint i_level,
    const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_last_ones) const {
    if (ptr_last_ones->size() != 3) {
        LogManager::LogError("size of ptr_last_ones should be 3 in MpiManager::GetNLevelCorrespondingOnes3D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    ptr_last_ones->at(kXIndex) =
        bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_)>>(kSFBitsetBit - i_level * 3);
    ptr_last_ones->at(kYIndex) =
        bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_)>>(kSFBitsetBit - i_level * 3);
    ptr_last_ones->at(kZIndex) =
        bitset_aux3d.k0SFBitsetTakeZRef_.at(bitset_aux3d.kRefCurrent_)>>(kSFBitsetBit - i_level * 3);
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[in] bitset_aux3d class manage space filling curves.
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 3
 */
void MpiManager::GetMinM1AtGivenLevel3D(const DefAmrIndexUint i_level,
    std::array<DefAmrIndexLUint, 3> indices_min, const SFBitsetAux3D& bitset_aux3d,
    std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const {
    if (ptr_min_m1_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 3 in MpiManager::GetMinM1AtGivenLevel3D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
     bitset_aux3d.SFBitsetEncoding({indices_min[kXIndex], 0, 0}));
    ptr_min_m1_bitsets->at(kXIndex) = bitset_aux3d.FindXNeg(bitset_tmp);
    bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
     bitset_aux3d.SFBitsetEncoding({0,  indices_min[kYIndex], 0}));
    ptr_min_m1_bitsets->at(kYIndex) = bitset_aux3d.FindYNeg(bitset_tmp);
    bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
     bitset_aux3d.SFBitsetEncoding({0, 0, indices_min[kZIndex]}));
    ptr_min_m1_bitsets->at(kZIndex) = bitset_aux3d.FindZNeg(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of maximum indices plus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[in] bitset_aux3d class manage space filling curves.
 * @param[out] ptr_max_p1_bitsets a pointer to maximum indices plus
 * @throws ErrorType if the size of max_p1_bitsets is not 3
 */
void MpiManager::GetMaxP1AtGivenLevel3D(const DefAmrIndexUint i_level,
    std::array<DefAmrIndexLUint, 3> indices_max, const SFBitsetAux3D& bitset_aux3d,
    std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const {
    if (ptr_max_p1_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 3 in MpiManager::GetMaxP1AtGivenLevel3D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
       bitset_aux3d.SFBitsetEncoding({indices_max[kXIndex], 0, 0}));
    ptr_max_p1_bitsets->at(kXIndex) = bitset_aux3d.FindXPos(bitset_tmp);
    bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
     bitset_aux3d.SFBitsetEncoding({0, indices_max[kYIndex], 0}));
    ptr_max_p1_bitsets->at(kYIndex) = bitset_aux3d.FindYPos(bitset_tmp);
    bitset_tmp = bitset_aux3d.SFBitsetToNHigherLevel(i_level,
     bitset_aux3d.SFBitsetEncoding({0, 0, indices_max[kZIndex]}));
    ptr_max_p1_bitsets->at(kZIndex) = bitset_aux3d.FindZPos(bitset_tmp);
}
/**
 * @brief function to check if a given node is on the interface of partitioned blocks
 * @param[in] i_level refinement level of input node.
 * @param[in] code_min minimum space fill code of current rank and specified refinement level.
 * @param[in] code_max maximum space fill code of current rank and specified refinement level.
 * @param[in] bitset_in space filling code of a given node.
 * @param[in] bitset_aux3d class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[in] bitset_level_ones bitsets of current refinement level excluding background space filling code.
 * @param[out] partitioned_interface_background  background nodes on the partitioned interface.
 */
bool MpiManager::CheckNodeOnOuterBoundaryOfBackgroundCell3D(DefAmrIndexUint i_level,
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const DefSFBitset bitset_in, const SFBitsetAux3D& bitset_aux3d,
    const std::vector<DefSFBitset>& domain_min_m1_n_level, const std::vector<DefSFBitset>& domain_max_p1_n_level,
    const std::vector<DefSFBitset>& bitset_level_ones,
    const DefMap<DefAmrIndexUint>& partitioned_interface_background) const {
    if ((bitset_in & bitset_level_ones.at(kXIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kXIndex)) == bitset_level_ones.at(kXIndex)
        || (bitset_in & bitset_level_ones.at(kYIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kYIndex)) == bitset_level_ones.at(kYIndex)
        || (bitset_in & bitset_level_ones.at(kZIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kZIndex)) == bitset_level_ones.at(kZIndex)) {
        DefSFBitset bitset_background = bitset_aux3d.SFBitsetToNLowerLevel(i_level, bitset_in), bitset_tmp;
        if (partitioned_interface_background.find(bitset_background) != partitioned_interface_background.end()) {
            std::array<DefSFBitset, 27> array_neighbors;
            bitset_aux3d.SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
            DefSFCodeToUint code;
            for (unsigned int i = 1; i < 27; ++i) {
                code = array_neighbors.at(i).to_ullong();
                if ((code < code_min || code > code_max)
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_min_m1_n_level.at(kXIndex))  // node is not x_min - 1
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_max_p1_n_level.at(kXIndex))
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_min_m1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_max_p1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeZRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_min_m1_n_level.at(kZIndex))
                    && ((array_neighbors.at(i) & bitset_aux3d.k0SFBitsetTakeZRef_.at(bitset_aux3d.kRefCurrent_))
                    != domain_max_p1_n_level.at(kZIndex))) {
                    return true;
                }
            }
        }
    }
    return false;
}
/**
 * @brief function to search for the ghost layers near a given node based on min and max space fill codes.
 * @param[in] sfbitset_in space fill code of the given node
 * @param[in] num_of_ghost_layers number of ghost layers
 * @param[in] code_min the minimum space fill codes of current rank and specified refinement level.
 * @param[in] code_max the maximum space fill codes of current rank and specified refinement level.
 * @param[in] flag_ini flag for initialization.
 * @param[in] bitset_aux3d class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[out] ptr_map_ghost_layer pointer to nodes on ghost layers near the given node.
 * @throws None
 */
void MpiManager::SearchForGhostLayerForMinNMax3D(const DefSFBitset sfbitset_in,
    const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const DefAmrIndexUint flag_ini, const SFBitsetAux3D& bitset_aux3d,
    const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level,
    DefMap<DefAmrIndexUint>* const ptr_map_ghost_layer) const {
    ptr_map_ghost_layer->clear();
    DefSFCodeToUint code_tmp;
    DefSFBitset sfbitset_tmp_y, sfbitset_tmp_x, sfbitset_tmp_z = sfbitset_in;
    // negative z direction
    for (DefAmrIndexUint iz = 0; iz <= num_of_ghost_layers; ++iz) {
        if ((sfbitset_tmp_z&bitset_aux3d.k0SFBitsetTakeZRef_.at(bitset_aux3d.kRefCurrent_))
         != domain_min_m1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefAmrIndexUint iy = 0; iy <= num_of_ghost_layers; ++iy) {
                if ((sfbitset_tmp_y&bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = bitset_aux3d.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = bitset_aux3d.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = bitset_aux3d.FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefAmrIndexUint iy = 0; iy < num_of_ghost_layers; ++iy) {
                sfbitset_tmp_y = bitset_aux3d.FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = bitset_aux3d.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = bitset_aux3d.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }
        } else {
            break;
        }
        sfbitset_tmp_z = bitset_aux3d.FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    for (DefAmrIndexUint iz = 0; iz < num_of_ghost_layers; ++iz) {
        sfbitset_tmp_z = bitset_aux3d.FindZPos(sfbitset_tmp_z);
        if ((sfbitset_tmp_z&bitset_aux3d.k0SFBitsetTakeZRef_.at(bitset_aux3d.kRefCurrent_))
         != domain_max_p1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefAmrIndexUint iy = 0; iy <= num_of_ghost_layers; ++iy) {
                if ((sfbitset_tmp_y&bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = bitset_aux3d.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = bitset_aux3d.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = bitset_aux3d.FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefAmrIndexUint iy = 0; iy < num_of_ghost_layers; ++iy) {
                sfbitset_tmp_y = bitset_aux3d.FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&bitset_aux3d.k0SFBitsetTakeYRef_.at(bitset_aux3d.kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = bitset_aux3d.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = bitset_aux3d.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&bitset_aux3d.k0SFBitsetTakeXRef_.at(bitset_aux3d.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_tmp_x.to_ullong();
                            if (code_tmp > code_max || code_tmp < code_min) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }
        } else {
            break;
        }
    }
}
/**
* @brief   function to calculate space filling code for mpi partition.
* @param[in] bitset_domain_min minimum space filling code of the computational domain
* @param[in] bitset_domain_max maximum space filling code of the computational domain
* @param[in] vec_cost computational cost from 0 to n refinement levels
* @param[in] vec_sfbitset space-filling codes for nodes of entire mesh.
* @param[in] bitset_aux3d class manage 3D space filling curves.
* @param[out] ptr_bitset_min pointer to minimum space filling code for each rank.
* @param[out] ptr_bitset_max pointer to maximum space filling code for each rank.
*/
void MpiManager::TraverseBackgroundForPartitionRank0(
    const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
    const std::vector<DefAmrIndexLUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
    const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_bitset_min,
    std::vector<DefSFBitset>* const ptr_bitset_max) const {
    DefMap<DefAmrUint> background_occupied;
    DefSFBitset bitset_background;
    DefAmrIndexUint max_level = DefAmrIndexUint(vec_cost.size()) - 1;
    DefAmrUint bk_cost =  vec_cost.at(0);
    std::array<DefAmrIndexLUint, 3> indices_min,  indices_max;
    bitset_aux3d.SFBitsetComputeIndices(bitset_domain_min, &indices_min);
    bitset_aux3d.SFBitsetComputeIndices(bitset_domain_max, &indices_max);
    DefAmrIndexLUint num_background_nodes = (indices_max[kXIndex] - indices_min[kXIndex] + 1)
        *(indices_max[kYIndex] - indices_min[kYIndex] + 1)*(indices_max[kZIndex] - indices_min[kZIndex] + 1);
    DefAmrIndexLUint sum_load = num_background_nodes * bk_cost;
    // add computational cost of refined nodes
    for (DefAmrIndexUint i_level = max_level; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        DefAmrUint node_cost = vec_cost.at(i_level);
        for (const auto& iter_low : vec_sfbitset.at(i_level)) {
            bitset_background = bitset_aux3d.SFBitsetToNLowerLevel(i_level_lower, iter_low.first);
            if (background_occupied.find(bitset_background) ==
                background_occupied.end()) {
                background_occupied.insert({ bitset_background, node_cost });
                sum_load += (node_cost - bk_cost);
            } else {
                background_occupied.at(bitset_background) += node_cost;
                sum_load += node_cost;
            }
        }
    }
    // calculate loads on each rank
    const int num_ranks = num_of_ranks_;
    DefAmrIndexLUint ave_load = static_cast<DefAmrIndexLUint>(sum_load / num_ranks) + 1;
    DefAmrIndexLUint load_rank0 = sum_load - (num_ranks - 1) * ave_load;
    std::vector<DefAmrIndexLUint> rank_load(num_ranks, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0, assuming lower than other ranks
    // traverse background nodes
    ptr_bitset_min->resize(num_ranks);
    ptr_bitset_max->resize(num_ranks);
    DefAmrIndexLUint load_count = 0;
    int status;
    int i_rank = 0;
    ptr_bitset_min->at(i_rank) = 0;
    ptr_bitset_max->back() = bitset_domain_max;
    DefAmrIndexLUint cost_background = vec_cost.at(0);
    std::array<DefAmrIndexLUint, 3> indices(indices_min);
    DefSFCodeToUint i_code = bitset_domain_min.to_ullong();
    DefSFBitset bitset_temp = static_cast<DefSFBitset>(i_code);
    for (DefAmrIndexLUint i_node = 0; i_node < num_background_nodes - 1; ++i_node) {
        if (load_count >= rank_load.at(i_rank)) {
            load_count = 0;
            ++i_rank;
            ptr_bitset_min->at(i_rank) = bitset_temp;
        }
        if (background_occupied.find(bitset_temp)
            == background_occupied.end()) {
            load_count += bk_cost;
        } else {
            load_count += background_occupied.at(bitset_temp);
        }
        if (load_count >= rank_load.at(i_rank)) {
            ptr_bitset_max->at(i_rank) = bitset_temp;
        }
        //  reset i_code if indices exceed domain
        ++i_code;
        status = bitset_aux3d.ResetIndicesExceedingDomain(indices_min, indices_max, &i_code, &bitset_temp);
#ifdef DEBUG_CHECK_GRID
        if (status) {
            LogManager::LogError("iterations exceed the maximum when space filling code exceed domain boundary "
                "in ResetIndicesExceedingDomain in MpiManager::TraverseBackgroundForPartition");
        }
#endif
    }
}
/**
 * @brief function to find interface of partitioned blocks
 * @param[in] bitset_min   minimum space filling code for each partition.
 * @param[in] bitset_max   maximum space filling code for each partition.
 * @param[in] array_domain_min minimum code of the computational domain.
 * @param[in] array_domain_max minimum code of the computational domain.
 * @param[in] bitset_aux3d class manage 3D space filling curves.
 * @param[out] ptr_partition_interface_background pointer to nodes on the interface of partitioned blocks at the background level. 
 */
void MpiManager::FindInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
    const DefSFBitset& bitset_max, const std::array<DefAmrIndexLUint, 3>& array_domain_min,
    const std::array<DefAmrIndexLUint, 3>& array_domain_max, const SFBitsetAux3D& bitset_aux3d,
    DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const {
    const DefSFCodeToUint code_min =  bitset_min.to_ullong(), code_max =  bitset_max.to_ullong();
    DefSFCodeToUint code_tmp = code_max + 1;
    DefAmrIndexUint block_level = 0;
    DefAmrIndexLUint block_length = 1 << block_level;  // block size of space filling code
    DefSFCodeToUint code_cri = ((code_tmp + 1) / 8)*block_length*block_length*block_length *8;

    DefSFCodeToUint code_max_criterion, code_remain;
    DefSFCodeToUint code_min_current;
    while (code_cri >= code_min) {
        code_max_criterion = code_tmp;
        code_min_current = code_max/ block_length/ block_length/ block_length;
        // though using FindPartitionRemainMax standalone gives the same result,
        // using FindPartitionBlocksMax additionally will take less iterations
        if (code_tmp - (code_tmp  - 1)%8 - 1 < code_min_current) {
            bitset_aux3d.FindPartitionRemainMax(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            bitset_aux3d.FindPartitionBlocksMax(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        }
        code_tmp/=8;
        block_level += 1;
        block_length = 1 << block_level;
        code_remain = code_tmp%8;
        if (code_remain == 0) {
            code_cri = (code_tmp - 8)/8;
        } else {
            code_cri = code_tmp /8;
        }
        code_cri = code_cri*block_length*block_length*block_length*8;
    }
    if (code_max_criterion%8 == 0) {
        code_tmp = code_max_criterion - 8;
    } else {
        code_tmp = code_max_criterion - code_max_criterion%8;
    }
    code_max_criterion = code_tmp * block_length * block_length * block_length/ 8;

    block_level = 0;
    block_length = 1;
    code_tmp = code_min;
    code_cri = code_min;
    DefSFCodeToUint code_max_current;
    while (code_cri < code_max_criterion) {
        code_max_current = code_max/ block_length/ block_length/ block_length;
        //  though using FindPartitionRemainMin standalone gives the same result,
        //  using FindPartitionBlocksMin additionally will take less iterations
        if (code_tmp + (8 - code_tmp%8) > code_max_current) {
            bitset_aux3d.FindPartitionRemainMin(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            bitset_aux3d.FindPartitionBlocksMin(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        }
        block_level += 1;
        block_length = 1 << block_level;
        code_remain = 8 - (code_tmp % 8);
        code_tmp = (code_tmp + code_remain) / 8;
        code_cri = code_tmp *block_length *block_length *block_length;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
