//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager2d_mpi.cpp
* @author Zhengliang Liu
* @date  2022-5-23
* @brief functions used to manage grid related processes when mpi is enabled.
*/
#include "../defs_libs.h"
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief 
 * @param bitset_min 
 * @param bitset_max 
 * @param ptr_partition_interface_background 
 * @note ptr_partition_interface_background include nodes outside the domain boundary,
 *       which should not be used in communication. These will be neglected during instantiate nodes.
 */
void GridManager2D::FindBlockInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
    const DefSFBitset& bitset_max, const std::vector<DefUint>& code_domain_min,
    const std::vector<DefUint>& code_domain_max, DefMap<DefUint>* const ptr_partition_interface_background) {
    const DefSFCodeToUint code_min =  bitset_min.to_ullong(), code_max =  bitset_max.to_ullong();
    DefSFCodeToUint code_tmp = code_max + 1;
    DefUint block_length = 1;  // block size of space filling code
    DefSFCodeToUint code_cri = ((code_tmp + 1) / 4) *block_length *block_length *4;

    std::array<DefUint, 2> array_domain_min = {code_domain_min.at(kXIndex), code_domain_min.at(kYIndex)},
        array_domain_max = {code_domain_max.at(kXIndex), code_domain_max.at(kYIndex)};
    DefSFCodeToUint code_max_criterion, code_remain;
    DefSFCodeToUint code_min_current, code_max_current;
    while (code_cri >= code_min) {
        code_max_criterion = code_tmp;
        code_min_current = code_max/ block_length/ block_length;
        // though using FindPartitionRemainMax standalone gives the same result,
        // using FindPartitionBlocksMax additionally will take less iterations
        if (code_tmp - (code_tmp  - 1)%4 - 1 < code_min_current) {
            FindPartitionRemainMax(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            FindPartitionBlocksMax(code_tmp, block_length, code_min, code_max,
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
            FindPartitionRemainMin(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            FindPartitionBlocksMin(code_tmp, block_length, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        }
        block_length *= 2;
        code_remain = 4 - (code_tmp % 4);
        code_tmp = (code_tmp + code_remain) / 4;
        code_cri = code_tmp *block_length *block_length;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
