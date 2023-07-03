//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager2d_mpi.cpp
* @author Zhengliang Liu
* @date  2022-5-23
* @brief functions used to manage grid related processes when mpi is enabled.
*/
#include "../defs_libs.h"
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
void GridManager3D::FindBlockInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
    const DefSFBitset& bitset_max, const std::vector<DefUint>& code_domain_min,
    const std::vector<DefUint>& code_domain_max, DefMap<DefUint>* const ptr_partition_interface_background) {
    const DefSFCodeToUint code_min =  bitset_min.to_ullong(), code_max =  bitset_max.to_ullong();
    DefSFCodeToUint code_tmp = code_max + 1;
    DefUint block_level = 0;
    DefUint block_length = 1 << block_level;  // block size of space filling code
    DefSFCodeToUint code_cri = ((code_tmp + 1) / 8)*block_length*block_length*block_length *8;

    std::array<DefUint, 3> array_domain_min = {code_domain_min.at(kXIndex),
        code_domain_min.at(kYIndex), code_domain_min.at(kZIndex)},
        array_domain_max = {code_domain_max.at(kXIndex), code_domain_max.at(kYIndex), code_domain_max.at(kZIndex)};
    DefSFCodeToUint code_max_criterion, code_remain;
    DefSFCodeToUint code_min_current;
    while (code_cri >= code_min) {
        code_max_criterion = code_tmp;
        code_min_current = code_max/ block_length/ block_length/ block_length;
        // though using FindPartitionRemainMax standalone gives the same result,
        // using FindPartitionBlocksMax additionally will take less iterations
        if (code_tmp - (code_tmp  - 1)%8 - 1 < code_min_current) {
            FindPartitionRemainMax(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        } else {
            FindPartitionBlocksMax(code_tmp, block_level, code_min, code_max,
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
    code_tmp = code_min;
    code_cri = code_min;
    DefSFCodeToUint code_max_current;
    while (code_cri < code_max_criterion) {
        code_max_current = code_max/ block_length/ block_length/ block_length;
        // though using FindPartitionRemainMin standalone gives the same result,
        // using FindPartitionBlocksMin additionally will take less iterations
        // if (code_tmp + (8 - code_tmp%8) > code_max_current) {
            FindPartitionRemainMin(code_tmp, block_level, code_min, code_max,
            array_domain_min, array_domain_max, ptr_partition_interface_background);
        // } else {
        //     FindPartitionBlocksMin(code_tmp, block_level, code_min, code_max,
        //     array_domain_min, array_domain_max, ptr_partition_interface_background);
        // }
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