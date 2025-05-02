//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
 * @brief function to check if a given node is on the interface of partitioned blocks
 * @param[in] i_level refinement level of input node.
 * @param[in] code_min minimum space fill code of current rank and specified refinement level.
 * @param[in] code_max maximum space fill code of current rank and specified refinement level.
 * @param[in] bitset_in space filling code of a given node.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[in] bitset_level_ones bitsets of current refinement level excluding background space filling code.
 * @param[out] partitioned_interface_background  background nodes on the partitioned interface.
 * @return 0 is not on the partitioned interface; 1 is on the lower interface; and 2 on the upper interface
 */
int MpiManager::CheckNodeOnPartitionInterface3D(DefInt i_level,
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const DefSFBitset bitset_in, const SFBitsetAux3D& sfbitset_aux,
    const std::vector<DefSFBitset>& domain_min_m1_n_level, const std::vector<DefSFBitset>& domain_max_p1_n_level,
    const std::vector<DefSFBitset>& bitset_level_ones,
    const DefMap<DefInt>& partitioned_interface_background) const {
    // noting that some neighbors of a nodes can be less than the minimum and some are greater than the maximum
    // thus |= operator other than a single return value is used to take this into consideration
    const std::array<DefSFBitset, 2>& take_xref = sfbitset_aux.GetTakeXRef(),
        take_yref =  sfbitset_aux.GetTakeYRef(), take_zref = sfbitset_aux.GetTakeZRef();
    int interface_status = 0;
    if ((bitset_in & bitset_level_ones.at(kXIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kYIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kZIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kXIndex)) == bitset_level_ones.at(kXIndex)
        || (bitset_in & bitset_level_ones.at(kYIndex)) == bitset_level_ones.at(kYIndex)
        || (bitset_in & bitset_level_ones.at(kZIndex)) == bitset_level_ones.at(kZIndex)) {
        DefSFBitset bitset_background = sfbitset_aux.SFBitsetToNLowerLevel(i_level, bitset_in), bitset_tmp;
        if (partitioned_interface_background.find(bitset_background) != partitioned_interface_background.end()) {
            std::array<DefSFBitset, 27> array_neighbors;
            sfbitset_aux.SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
            DefSFCodeToUint code;
            for (unsigned int i = 1; i < 27; ++i) {
                code = sfbitset_aux.SFBitsetoSFCode(array_neighbors.at(i));
                if (code < code_min
                    && ((array_neighbors.at(i) & take_xref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kXIndex))  // node is not x_min - 1
                    && ((array_neighbors.at(i) & take_xref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kXIndex))
                    && ((array_neighbors.at(i) & take_yref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & take_yref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & take_zref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kZIndex))
                    && ((array_neighbors.at(i) & take_zref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kZIndex))) {
                    interface_status |= 1;
                } else if (code > code_max
                    && ((array_neighbors.at(i) & take_xref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kXIndex))  // node is not x_min - 1
                    && ((array_neighbors.at(i) & take_xref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kXIndex))
                    && ((array_neighbors.at(i) & take_yref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & take_yref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & take_zref.at(sfbitset_aux.kRefCurrent_))
                    != domain_min_m1_n_level.at(kZIndex))
                    && ((array_neighbors.at(i) & take_zref.at(sfbitset_aux.kRefCurrent_))
                    != domain_max_p1_n_level.at(kZIndex))) {
                    interface_status |= 2;
                }
            }
        }
    }
    return interface_status;
}
/**
 * @brief function to search for the ghost layers near a given node based on min and max space fill codes.
 * @param[in] sfbitset_in space fill code of the given node.
 * @param[in] num_of_ghost_layers number of ghost layers
 * @param[in] code_bound the minimum and maximum space fill code.
 * @param[in] ptr_func_compare pointer to function to check if code is less than or greater than the bounds.
 * @param[in] flag_ini flag for initialization.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[out] ptr_map_ghost_layer pointer to nodes on ghost layers near the given node.
 */
void MpiManager::SearchForGhostLayerForMinNMax3D(const DefSFBitset sfbitset_in,
    const DefInt num_of_ghost_layers, const DefSFCodeToUint code_bound,
    bool (MpiManager::*ptr_func_compare)(const DefSFCodeToUint, const DefSFCodeToUint) const,
    const DefInt flag_ini, const SFBitsetAux3D& sfbitset_aux,
    const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level,
    DefMap<DefInt>* const ptr_map_ghost_layer) const {
    DefSFCodeToUint code_tmp;
    DefSFBitset sfbitset_tmp_y, sfbitset_tmp_x, sfbitset_tmp_z = sfbitset_in;
    const std::array<DefSFBitset, 2>& take_xref = sfbitset_aux.GetTakeXRef(),
        take_yref =  sfbitset_aux.GetTakeYRef(), take_zref = sfbitset_aux.GetTakeZRef();
    // negative z direction
    for (DefInt iz = 0; iz <= num_of_ghost_layers; ++iz) {
        if ((sfbitset_tmp_z&take_zref.at(sfbitset_aux.kRefCurrent_))
         != domain_min_m1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefInt iy = 0; iy <= num_of_ghost_layers; ++iy) {
                if ((sfbitset_tmp_y&take_yref.at(sfbitset_aux.kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = sfbitset_aux.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = sfbitset_aux.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = sfbitset_aux.FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefInt iy = 0; iy < num_of_ghost_layers; ++iy) {
                sfbitset_tmp_y = sfbitset_aux.FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&take_yref.at(sfbitset_aux.kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = sfbitset_aux.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = sfbitset_aux.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
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
        sfbitset_tmp_z = sfbitset_aux.FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    for (DefInt iz = 0; iz < num_of_ghost_layers; ++iz) {
        sfbitset_tmp_z = sfbitset_aux.FindZPos(sfbitset_tmp_z);
        if ((sfbitset_tmp_z&take_zref.at(sfbitset_aux.kRefCurrent_))
         != domain_max_p1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefInt iy = 0; iy <= num_of_ghost_layers; ++iy) {
                if ((sfbitset_tmp_y&take_yref.at(sfbitset_aux.kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = sfbitset_aux.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = sfbitset_aux.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = sfbitset_aux.FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefInt iy = 0; iy < num_of_ghost_layers; ++iy) {
                sfbitset_tmp_y = sfbitset_aux.FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&take_yref.at(sfbitset_aux.kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= num_of_ghost_layers; ++ix) {
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
                                ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                            }
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = sfbitset_aux.FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < num_of_ghost_layers; ++ix) {
                        sfbitset_tmp_x = sfbitset_aux.FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&take_xref.at(sfbitset_aux.kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            code_tmp = sfbitset_aux.SFBitsetoSFCode(sfbitset_tmp_x);
                            if ((this->*ptr_func_compare)(code_tmp, code_bound)) {
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
* @param[in] sfbitset_aux class manage 3D space filling curves.
* @param[out] ptr_bitset_min pointer to minimum space filling code for each rank.
* @param[out] ptr_bitset_max pointer to maximum space filling code for each rank.
*/
void MpiManager::IniTraverseBackgroundForPartitionRank0(
    const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
    const std::vector<DefInt>& vec_cost, const std::vector<DefMap<DefInt>>& vec_sfbitset,
    const SFBitsetAux3D& sfbitset_aux, std::vector<DefSFBitset>* const ptr_bitset_min,
    std::vector<DefSFBitset>* const ptr_bitset_max) const {
    DefMap<DefInt> background_occupied;
    DefSFBitset bitset_background;
    DefInt max_level = DefInt(vec_cost.size()) - 1;
    DefInt bk_cost =  vec_cost.at(0);
    std::array<DefAmrLUint, 3> indices_min,  indices_max;
    sfbitset_aux.SFBitsetComputeIndices(bitset_domain_min, &indices_min);
    sfbitset_aux.SFBitsetComputeIndices(bitset_domain_max, &indices_max);
    DefAmrLUint num_background_nodes = (indices_max[kXIndex] - indices_min[kXIndex] + 1)
        *(indices_max[kYIndex] - indices_min[kYIndex] + 1)*(indices_max[kZIndex] - indices_min[kZIndex] + 1);
    DefAmrLUint sum_load = num_background_nodes * bk_cost;
    // add computational cost of refined nodes
    for (DefInt i_level = max_level; i_level > 0; --i_level) {
        DefInt i_level_lower = i_level - 1;
        DefInt node_cost = vec_cost.at(i_level);
        for (const auto& iter_low : vec_sfbitset.at(i_level)) {
            bitset_background = sfbitset_aux.SFBitsetToNLowerLevel(i_level_lower, iter_low.first);
            if (sum_load + node_cost > (std::numeric_limits<DefAmrLUint>::max)()) {
                LogManager::LogError("computational load exceeds the maximum of DefAmrLUint in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
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
    DefAmrLUint ave_load = static_cast<DefAmrLUint>(sum_load / num_ranks) + 1;
    DefAmrLUint load_rank0 = sum_load - (num_ranks - 1) * ave_load;
    std::vector<DefAmrLUint> rank_load(num_ranks, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0, assuming lower than other ranks
    // traverse background nodes
    ptr_bitset_min->resize(num_ranks);
    ptr_bitset_max->resize(num_ranks);
    DefAmrLUint load_count = 0;
    int status;
    int i_rank = 0;
    ptr_bitset_min->at(i_rank) = 0;
    ptr_bitset_max->back() = bitset_domain_max;
    std::array<DefAmrLUint, 3> indices(indices_min);
    DefSFCodeToUint i_code = sfbitset_aux.SFBitsetoSFCode(bitset_domain_min);
    DefSFBitset sfbitset_tmp = static_cast<DefSFBitset>(i_code);
    for (DefAmrLUint i_node = 0; i_node < num_background_nodes - 1; ++i_node) {
        if (load_count >= rank_load.at(i_rank)) {
            load_count = 0;
            ++i_rank;
            ptr_bitset_min->at(i_rank) = sfbitset_tmp;
        }
        if (background_occupied.find(sfbitset_tmp)
            == background_occupied.end()) {
            load_count += bk_cost;
        } else {
            load_count += background_occupied.at(sfbitset_tmp);
        }
        if (load_count >= rank_load.at(i_rank)) {
            ptr_bitset_max->at(i_rank) = sfbitset_tmp;
        }
        //  reset i_code if indices exceed domain
        ++i_code;
        status = sfbitset_aux.ResetIndicesExceedingDomain(indices_min, indices_max, &i_code, &sfbitset_tmp);
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
 * @param[in] code_min   minimum space filling code for each partition.
 * @param[in] code_max   maximum space filling code for each partition.
 * @param[in] array_domain_min minimum code of the computational domain.
 * @param[in] array_domain_max minimum code of the computational domain.
 * @param[in] sfbitset_aux class manage 3D space filling curves.
 * @param[out] ptr_partition_interface_background pointer to nodes on the interface of partitioned blocks at the background level. 
 */
void MpiManager::IniFindInterfaceForPartitionFromMinNMax(const DefSFCodeToUint& code_min,
    const DefSFCodeToUint& code_max, const std::array<DefAmrLUint, 3>& array_domain_min,
    const std::array<DefAmrLUint, 3>& array_domain_max, const SFBitsetAux3D& sfbitset_aux,
    DefMap<DefInt>* const ptr_partition_interface_background) const {
    DefSFCodeToUint code_tmp = code_max + 1;
    DefInt block_level = 0;
    DefAmrLUint block_length = 1 << block_level;  // block size of space filling code
    DefSFCodeToUint code_cri = ((code_tmp + 1) / 8)*block_length*block_length*block_length *8;

    DefSFCodeToUint code_max_criterion, code_remain;
    DefSFCodeToUint code_min_current;
    while (code_cri >= code_min && code_cri > 0) {
        code_max_criterion = code_tmp;
        code_min_current = code_max/ block_length/ block_length/ block_length;
        // though using FindPartitionRemainMax standalone gives the same result,
        // using FindPartitionBlocksMax additionally will take less iterations
        if (code_tmp - (code_tmp  - 1)%8 - 1 < code_min_current) {
            sfbitset_aux.FindPartitionRemainMax(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            sfbitset_aux.FindPartitionBlocksMax(code_tmp, block_level, code_min, code_max,
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
            sfbitset_aux.FindPartitionRemainMin(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            sfbitset_aux.FindPartitionBlocksMin(code_tmp, block_level, code_min, code_max,
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
