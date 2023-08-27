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
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
namespace rootproject {
namespace amrproject {
/**
 * @brief function to get spacing fill codes whose bits are 1 for the given level and 0 for the background
 * @param[in] i_level given refinement level
 * @param[in] bitset_aux2d class manage space filling curves.
 * @param[out] ptr_last_ones pointer to spacing fill codes
 * @throws ErrorType if the size of last_ones is not 2
 */
void MpiManager::GetNLevelCorrespondingOnes2D(const DefAmrIndexUint i_level,
    const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_last_ones) const {
    if (ptr_last_ones->size() != 2) {
        LogManager::LogError("size of ptr_last_ones should be 2 in MpiManager::GetNLevelCorrespondingOnes2D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    ptr_last_ones->at(kXIndex) =
        bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_)>>(kSFBitsetBit - i_level * 2);
    ptr_last_ones->at(kYIndex) =
        bitset_aux2d.k0SFBitsetTakeYRef_.at(bitset_aux2d.kRefCurrent_)>>(kSFBitsetBit - i_level * 2);
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[in] bitset_aux2d class manage space filling curves.
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 2
 */
void MpiManager::GetMinM1AtGivenLevel2D(const DefAmrIndexUint i_level,
    std::array<DefAmrIndexLUint, 2> indices_min, const SFBitsetAux2D& bitset_aux2d,
    std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const {
    if (ptr_min_m1_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 2 in MpiManager::GetMinM1AtGivenLevel2D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = bitset_aux2d.SFBitsetToNHigherLevel(i_level,
     bitset_aux2d.SFBitsetEncoding({indices_min[kXIndex], 0}));
    ptr_min_m1_bitsets->at(kXIndex) = bitset_aux2d.FindXNeg(bitset_tmp);
    bitset_tmp = bitset_aux2d.SFBitsetToNHigherLevel(i_level, bitset_aux2d.SFBitsetEncoding({0, indices_min[kYIndex]}));
    ptr_min_m1_bitsets->at(kYIndex) = bitset_aux2d.FindYNeg(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of maximum indices plus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[in] bitset_aux2d class manage space filling curves.
 * @param[out] ptr_max_p1_bitsets a pointer to maximum indices plus
 * @throws ErrorType if the size of max_p1_bitsets is not 2
 */
void MpiManager::GetMaxP1AtGivenLevel2D(const DefAmrIndexUint i_level,
    std::array<DefAmrIndexLUint, 2> indices_max, const SFBitsetAux2D& bitset_aux2d,
    std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const {
    if (ptr_max_p1_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 2 in MpiManager::GetMaxP1AtGivenLevel2D in "
         + std::string(__FILE__) + "at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = bitset_aux2d.SFBitsetToNHigherLevel(i_level,
     bitset_aux2d.SFBitsetEncoding({indices_max[kXIndex], 0}));
    ptr_max_p1_bitsets->at(kXIndex) = bitset_aux2d.FindXPos(bitset_tmp);
    bitset_tmp = bitset_aux2d.SFBitsetToNHigherLevel(i_level,
     bitset_aux2d.SFBitsetEncoding({0, indices_max[kYIndex]}));
    ptr_max_p1_bitsets->at(kYIndex) = bitset_aux2d.FindYPos(bitset_tmp);
}
/**
 * @brief function to check if a given node is on the interface of partitioned blocks
 * @param[in] i_level refinement level of input node.
 * @param[in] code_min minimum space fill code of current rank and specified refinement level.
 * @param[in] code_max maximum space fill code of current rank and specified refinement level.
 * @param[in] bitset_in space filling code of a given node.
 * @param[in] bitset_aux2d class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[in] bitset_level_ones bitsets of current refinement level excluding background space filling code.
 * @param[out] partitioned_interface_background  background nodes on the partitioned interface.
 */
bool MpiManager::CheckNodeOnOuterBoundaryOfBackgroundCell2D(DefAmrIndexUint i_level,
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const DefSFBitset bitset_in, const SFBitsetAux2D& bitset_aux2d,
    const std::vector<DefSFBitset>& domain_min_m1_n_level, const std::vector<DefSFBitset>& domain_max_p1_n_level,
    const std::vector<DefSFBitset>& bitset_level_ones,
    const DefMap<DefAmrIndexUint>& partitioned_interface_background) const {
    // At least one index is the minimum or maximum of the background cell,
    // or the node is inside the background cell, not on the edge.
    if ((bitset_in & bitset_level_ones.at(kXIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kXIndex)) == bitset_level_ones.at(kXIndex)
        || (bitset_in & bitset_level_ones.at(kYIndex)) == 0
        || (bitset_in & bitset_level_ones.at(kYIndex)) == bitset_level_ones.at(kYIndex)) {
        DefSFBitset bitset_background = bitset_aux2d.SFBitsetToNLowerLevel(i_level, bitset_in), bitset_tmp;
        if (partitioned_interface_background.find(bitset_background) != partitioned_interface_background.end()) {
            std::array<DefSFBitset, 9> array_neighbors;
            bitset_aux2d.SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
            DefSFCodeToUint code;
            for (unsigned int i = 1; i < 9; ++i) {
                code = array_neighbors.at(i).to_ullong();
                if ((code < code_min || code > code_max)
                    && ((array_neighbors.at(i) & bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
                    != domain_min_m1_n_level.at(kXIndex))
                    && ((array_neighbors.at(i) & bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
                    != domain_max_p1_n_level.at(kXIndex))
                    && ((array_neighbors.at(i) & bitset_aux2d.k0SFBitsetTakeYRef_.at(bitset_aux2d.kRefCurrent_))
                    != domain_min_m1_n_level.at(kYIndex))
                    && ((array_neighbors.at(i) & bitset_aux2d.k0SFBitsetTakeYRef_.at(bitset_aux2d.kRefCurrent_))
                    != domain_max_p1_n_level.at(kYIndex))) {
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
 * @param[in] bitset_aux2d class manage space filling curves.
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[out] ptr_map_ghost_layer pointer to nodes on ghost layers near the given node.
 * @throws None
 */
void MpiManager::SearchForGhostLayerForMinNMax2D(const DefSFBitset sfbitset_in,
    const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const DefAmrIndexUint flag_ini, const SFBitsetAux2D& bitset_aux2d,
    const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level,
    DefMap<DefAmrIndexUint>* const ptr_map_ghost_layer) const {
    ptr_map_ghost_layer->clear();
    DefSFCodeToUint code_tmp;
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    // negative y direction
    for (DefAmrIndexUint iy = 0; iy <= num_of_ghost_layers; ++iy) {
        if ((sfbitset_tmp_y&bitset_aux2d.k0SFBitsetTakeYRef_.at(bitset_aux2d.kRefCurrent_))
         != domain_min_m1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                if ((sfbitset_tmp_x&bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    code_tmp = sfbitset_tmp_x.to_ullong();
                    if (code_tmp > code_max || code_tmp < code_min) {
                        ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                    }
                } else {
                    break;
                }
                sfbitset_tmp_x = bitset_aux2d.FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                sfbitset_tmp_x = bitset_aux2d.FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
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
        sfbitset_tmp_y = bitset_aux2d.FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    for (DefAmrIndexUint iy = 0; iy < num_of_ghost_layers; ++iy) {
        sfbitset_tmp_y = bitset_aux2d.FindYPos(sfbitset_tmp_y);
        if ((sfbitset_tmp_y&bitset_aux2d.k0SFBitsetTakeYRef_.at(bitset_aux2d.kRefCurrent_))
         != domain_max_p1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix <= num_of_ghost_layers; ++ix) {
                if ((sfbitset_tmp_x&bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    code_tmp = sfbitset_tmp_x.to_ullong();
                    if (code_tmp > code_max || code_tmp < code_min) {
                        ptr_map_ghost_layer->insert({sfbitset_tmp_x, flag_ini});
                    }
                } else {
                    break;
                }
                sfbitset_tmp_x = bitset_aux2d.FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix < num_of_ghost_layers; ++ix) {
                sfbitset_tmp_x = bitset_aux2d.FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&bitset_aux2d.k0SFBitsetTakeXRef_.at(bitset_aux2d.kRefCurrent_))
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
}
/**
* @brief   function to calculate space filling code for mpi partition.
* @param[in] bitset_domain_min minimum space filling code of the computational domain
* @param[in] bitset_domain_max maximum space filling code of the computational domain
* @param[in] vec_cost computational cost from 0 to n refinement levels
* @param[in] vec_sfbitset space-filling codes for nodes of entire mesh.
* @param[in] bitset_aux2d class manage 2D space filling curves.
* @param[out] ptr_bitset_min pointer to minimum space filling code for each rank.
* @param[out] ptr_bitset_max pointer to maximum space filling code for each rank.
*/
void MpiManager::TraverseBackgroundForPartitionRank0(
    const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
    const std::vector<DefAmrIndexLUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
    const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_bitset_min,
    std::vector<DefSFBitset>* const ptr_bitset_max) const {
    DefMap<DefAmrUint> background_occupied;
    DefSFBitset bitset_background;
    DefAmrIndexUint max_level = DefAmrIndexUint(vec_cost.size()) - 1;
    DefAmrUint bk_cost =  vec_cost.at(0);
    std::array<DefAmrIndexLUint, 2> indices_min,  indices_max;
    bitset_aux2d.SFBitsetComputeIndices(bitset_domain_min, &indices_min);
    bitset_aux2d.SFBitsetComputeIndices(bitset_domain_max, &indices_max);
    DefAmrIndexLUint num_background_nodes = (indices_max[kXIndex] - indices_min[kXIndex] + 1)
        *(indices_max[kYIndex] - indices_min[kYIndex] + 1);
    DefAmrIndexLUint sum_load = num_background_nodes * bk_cost;
    // add computational cost of refined nodes
    for (DefAmrIndexUint i_level = max_level; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        DefAmrUint node_cost = vec_cost.at(i_level);
        for (const auto& iter_low : vec_sfbitset.at(i_level)) {
            bitset_background = bitset_aux2d.SFBitsetToNLowerLevel(i_level_lower, iter_low.first);
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
    int num_ranks = num_of_ranks_;
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
            load_count += cost_background;
        } else {
            load_count += background_occupied.at(bitset_temp);
        }
        if (load_count >= rank_load.at(i_rank)) {
            ptr_bitset_max->at(i_rank) = bitset_temp;
        }
        //  reset i_code if indices exceed domain
        ++i_code;
        bitset_temp = static_cast<DefSFBitset>(i_code);
        status = bitset_aux2d.ResetIndicesExceedingDomain(indices_min, indices_max, &i_code, &bitset_temp);
#ifdef DEBUG_CHECK_GRID
        if (status) {
            LogManager::LogError("iterations exceed the maximum when space filling code exceed domain boundary"
                " in ResetIndicesExceedingDomain in MpiManager::TraverseBackgroundForPartition");
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
 * @param[in] bitset_aux2d class manage 2D space filling curves.
 * @param[out] ptr_partition_interface_background pointer to nodes on the interface of partitioned blocks at the background level. 
 */
void MpiManager::FindInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
    const DefSFBitset& bitset_max, const std::array<DefAmrIndexLUint, 2>& array_domain_min,
    const std::array<DefAmrIndexLUint, 2>& array_domain_max, const SFBitsetAux2D& bitset_aux2d,
    DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const {
    const DefSFCodeToUint code_min =  bitset_min.to_ullong(), code_max =  bitset_max.to_ullong();
    DefSFCodeToUint code_tmp = code_max + 1;
    DefAmrIndexLUint block_length = 1;  // block size of space filling code
    DefSFCodeToUint code_cri = ((code_tmp + 1) / 4) *block_length *block_length *4;

    DefSFCodeToUint code_max_criterion, code_remain;
    DefSFCodeToUint code_min_current, code_max_current;
    while (code_cri >= code_min && code_cri > 0) {
        code_max_criterion = code_tmp;
        code_min_current = code_max/ block_length/ block_length;
        // though using FindPartitionRemainMax standalone gives the same result,
        // using FindPartitionBlocksMax additionally will take less iterations
        if (code_tmp - (code_tmp  - 1)%4 - 1 < code_min_current) {
            bitset_aux2d.FindPartitionRemainMax(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            bitset_aux2d.FindPartitionBlocksMax(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        }
        code_tmp/=4;
        block_length *= 2;
        code_remain = code_tmp%4;
        if (code_remain == 0) {
            code_cri = (code_tmp - 4)/4;
        } else {
            code_cri = code_tmp /4;
        }
        code_cri = code_cri*block_length *block_length *4;
    }
    if (code_max_criterion%4 == 0) {
        code_tmp = code_max_criterion - 4;
    } else {
        code_tmp = code_max_criterion - code_max_criterion%4;
    }
    code_max_criterion = code_tmp * block_length * block_length / 4;

    block_length = 1;
    code_tmp = code_min;
    code_cri = code_min;
    while (code_cri < code_max_criterion) {
        code_max_current = code_max/ block_length/ block_length;
        // though using FindPartitionRemainMin standalone gives the same result,
        // using FindPartitionBlocksMin additionally will take less iterations
        if (code_tmp + (4 - code_tmp%4) > code_max_current) {
            bitset_aux2d.FindPartitionRemainMin(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            bitset_aux2d.FindPartitionBlocksMin(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        }
        block_length *= 2;
        code_remain = 4 - (code_tmp % 4);
        code_tmp = (code_tmp + code_remain) / 4;
        code_cri = code_tmp *block_length *block_length;
    }
}
void MpiManager::SearchPartitionInterfaceLayer2D(const DefAmrIndexUint flag0, const GridNode& grid_not_instance,
    const SFBitsetAux2D sfbitset_aux2d, const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info,
    const std::vector<DefMap<DefAmrIndexUint>> sfbitset_partition_interface,
    GridInfoInterface* const ptr_background_grid) {
    DefSizet max_level = sfbitset_partition_interface.size();
    DefSFBitset sfbitset_current, sfbitset_x, sfbitset_y;
    for (auto i_level = 1; i_level < max_level; ++i_level) {
        for (const auto iter_node : sfbitset_partition_interface.at(i_level)) {
            sfbitset_current = sfbitset_aux2d.SFBitsetToNHigherLevel(1, iter_node.first);
            mpi_communication_ghost_layers_.at(i_level).at(0).insert({sfbitset_current, flag0});
            sfbitset_x = sfbitset_aux2d.FindXPos(iter_node.first);
            if (sfbitset_partition_interface.at(i_level).find(sfbitset_x)
             != sfbitset_partition_interface.at(i_level).end()) {
                mpi_communication_ghost_layers_.at(i_level).at(0).insert(
                    {sfbitset_aux2d.FindXPos(sfbitset_current), flag0});
            }
            sfbitset_y = sfbitset_aux2d.FindXPos(iter_node.first);
            if (sfbitset_partition_interface.at(i_level).find(sfbitset_y)
             != sfbitset_partition_interface.at(i_level).end()) {
                mpi_communication_ghost_layers_.at(i_level).at(0).insert(
                    {sfbitset_aux2d.FindXPos(sfbitset_current), flag0});
            }
        }
    }

}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#endif  // ENABLE_MPI
