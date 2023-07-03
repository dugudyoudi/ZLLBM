//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.cpp
* @author Zhengliang Liu
* @brief functions used to manipulate morton codes.
* @date  2022-5-16
* @note  functions from other managers will be called.
*/
#include <string>
#include <array>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include "auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
 * @brief function to find interface of partitioned blocks by several blocks at given level based on the maximum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length - 1;
    code_remain = (code_in - 1) % 4;
    std::array<DefLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 0:
            for (DefLUint ix = 0; ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0; ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 1:
            for (DefLUint ix = 0; ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0; ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    &&CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 2:
            for (DefLUint ix = 0; ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0; ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1 + 1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1 + 1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 3:
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks by several blocks at given level based on the minimum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length;
    code_remain = code_in % 4;
    std::array<DefLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 3:
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 2:
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 1:
            for (DefLUint iy = 0;  iy < i_max1 + 1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1 + 1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 0:
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks one by one block at given level based on the maximum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    code_current = (code_in - 1) * block_length * block_length;
    code_remain = (code_in - 1) % 4;
    DefSFCodeToUint code_min_current = code_partition_min/ block_length/ block_length;
    DefSFBitset bitset_tmp = static_cast<DefSFBitset>(code_current);
    switch (code_remain) {
        case 0:
            if (code_in - 1 >= code_min_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 1:
            if (code_in - 1 >= code_min_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 2 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 2) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 2:
            if (code_in - 1 >= code_min_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 2 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 2) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 3 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 3) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 3:
            if (code_in - 1 >= code_min_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 2 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 2) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 3 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 3) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in - 4 >= code_min_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in - 4) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks one by one block at given level based on the minimum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length;
    code_remain = code_in % 4;
    std::array<DefLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    DefSFCodeToUint code_max_current = code_partition_max/ block_length/ block_length;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 3:
            if (code_in + 1 <= code_max_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 2:
            if (code_in + 1 <= code_max_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 2 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 1) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 1:
            if (code_in + 1 <= code_max_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 2 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 1) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 3 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 2) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        case 0:
            if (code_in + 1 <= code_max_current) {
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 2 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 1) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 3 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 2) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            if (code_in + 4 <= code_max_current) {
                bitset_tmp = static_cast<DefSFBitset>((code_in + 3) * block_length * block_length);
                FindPartitionOneBlock(bitset_tmp, block_length, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
            }
            break;
        default:
        break;
    }
}
/**
 * @brief function to find interface of partitioned blocks by on one block 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
         std::array<DefLUint, 2> indices;
    DefSFBitset bitset_tmp = bitset_corner;
    SFBitsetComputeIndices(bitset_corner, &indices);
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
        && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
        && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
            code_domain_min, code_domain_max)
        && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
        ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        DefLUint i_max1 = block_length - 1;
        for (DefLUint ix = 0;  ix < i_max1; ++ix) {
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            bitset_tmp = FindXPos(bitset_tmp);
            ++index_x;
        }
        for (DefLUint iy = 0;  iy < i_max1; ++iy) {
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            bitset_tmp = FindYPos(bitset_tmp);
            ++index_y;
        }
        for (DefLUint ix = 0;  ix < i_max1; ++ix) {
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            bitset_tmp = FindXNeg(bitset_tmp);
            --index_x;
        }
        for (DefLUint iy = 0;  iy < i_max1; ++iy) {
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            bitset_tmp = FindYNeg(bitset_tmp);
            --index_y;
        }
    }
}
/**
 * @brief function to check if a node is on the interface of partitioned blocks at the current rank
 * @param[in] sfbitset_in space filling code of input node.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @return true if the input node is on the interface of partitioned blocks at the current rank
 */
bool SFBitsetAux2D::CheckPartitionLimits(const DefSFBitset& sfbitset_in,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 2>& code_domain_min, const std::array<DefLUint, 2>& code_domain_max) const {
    std::array<DefSFBitset, 9> neighbours;
    DefSFCodeToUint code_neighbour;
    std::array<DefLUint, 2> indices;
    SFBitsetFindAllNeighbours(sfbitset_in, &neighbours);
    for (DefUint i = 1; i < 9; ++i) {
        code_neighbour = neighbours[i].to_ullong();
        if (code_neighbour < code_partition_min || code_neighbour > code_partition_max) {
            SFBitsetComputeIndices(neighbours[i], &indices);
            if (indices[kXIndex] < code_domain_min[kXIndex] || indices[kXIndex] > code_domain_max[kXIndex]
                || indices[kYIndex] < code_domain_min[kYIndex] || indices[kYIndex] > code_domain_max[kYIndex]) {
            } else {
                return true;
            }
        }
    }
    return false;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
void SFBitsetAux3D::FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_tmp;
    DefSFCodeToUint code_current;
    const DefUint block_length = 1 << block_level;
    DefSFCodeToUint block_length_cubic = block_length * block_length * block_length;
    code_remain = (code_in - 1) % 8;
    code_current = (code_in - 1 - code_remain) * block_length_cubic - code_remain;
    bitset_tmp = static_cast<DefSFBitset>(code_current);
    switch (code_remain) {
        case 0:
            FindPartitionOneBlock(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 1:
            FindPartition2BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 2:
            FindPartition3BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 3:
            FindPartition4BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 4:
            FindPartition5BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 5:
            FindPartition6BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 6:
            FindPartition7BlocksMax(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 7:
            FindPartitionOneBlock(bitset_tmp, block_level + 1, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks by several blocks at given level based on the minimum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    const DefUint block_length = 1 << block_level;
    DefSFCodeToUint block_length_cubic = block_length * block_length * block_length;
    code_remain = code_in % 4;
    DefSFCodeToUint code_current = code_in*block_length_cubic;
    std::array<DefLUint, 3> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 3:
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                      code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 2:
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 1:
            for (DefLUint iy = 0;  iy < i_max1 + 1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1 + 1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                   code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        case 0:
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                       code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
            if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max)
                && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                ptr_map_partitioned->insert({bitset_tmp, 0});
            }
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks one by one block at given level based on the maximum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    const DefUint block_length = 1 << block_level;
    const DefUint block_length_cubic =  block_length * block_length * block_length;
    code_current = (code_in - 1) * block_length_cubic;
    code_remain = (code_in - 1) % 8;
    DefSFCodeToUint code_min_current = code_partition_min/ block_length_cubic;
    DefSFBitset bitset_tmp = static_cast<DefSFBitset>(code_current);

    for (auto i = 0; i <= code_remain; ++i) {
        if (code_in - 1 - i>= code_min_current) {
            bitset_tmp = static_cast<DefSFBitset>((code_in - 1 - i) * block_length_cubic);
            FindPartitionOneBlock(bitset_tmp, block_level, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
        } else {
            break;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks one by one block at given level based on the minimum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    const DefUint block_length = 1 << block_level;
    const DefUint block_length_cubic =  block_length * block_length * block_length;
    code_current = code_in * block_length_cubic;
    code_remain = code_in % 8;
    DefSFCodeToUint code_max_current = code_partition_max/ block_length_cubic;
    DefSFBitset bitset_tmp = static_cast<DefSFBitset>(code_current);
    code_remain = 8 - code_remain;
    for (auto i = 1; i <= code_remain; ++i) {
        if (code_in + i <= code_max_current) {
            bitset_tmp = static_cast<DefSFBitset>((code_in - 1 + i) * block_length_cubic);
            FindPartitionOneBlock(bitset_tmp, block_level, code_partition_min, code_partition_max,
                    code_domain_min, code_domain_max, ptr_map_partitioned);
        } else {
            break;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by one block 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1;
    SFBitsetComputeIndices(bitset_corner, &indices);
    DefSFBitset bitset_tmp = bitset_corner, bitset_tmp2 = bitset_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by two blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition2BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by three blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition3BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        for (DefLUint iy = 0;  iy < i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0; ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0; ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy <= i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0; ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        for (DefLUint iy = 0;  iy < i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by four blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition4BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by five blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition5BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
         if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // lower vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // middle top surface
        std::array<DefLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindXPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices2[kXIndex];
            ++index_y;
        }
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices2[kXIndex];
            ++index_y;
        }
        // upper vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by six blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition6BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // lower vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // middle top surface
        std::array<DefLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices2[kXIndex];
            ++index_y;
        }
         // upper vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by seven blocks bounded by the maximum space fill code. 
 * @param[in] bitset_corner space filling code of left right conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition7BlocksMax(const DefSFBitset& bitset_min_corner, const DefUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max,
    DefMap<DefUint>* const ptr_map_partitioned) const {
    std::array<DefLUint, 3> indices;
    const DefUint block_length = 1 << block_level;
    DefLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp2 = FindXPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp2, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp2) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp2, 0});
        }
        bitset_tmp2 = FindYPos(bitset_tmp);
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp2, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp2) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp2, 0});
        }
    } else {
        // bottom surface
        for (DefLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        // lower vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0;  ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0;  iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // middle top surface
        std::array<DefLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindXPos(bitset_tmp);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices2[kXIndex];
            ++index_y;
        }
        // upper vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefLUint ix = 0; ix < i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            for (DefLUint iy = 0; iy < i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0; ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy <= i_max1; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYPos(bitset_tmp);
                ++index_y;
            }
            for (DefLUint ix = 0; ix < i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXNeg(bitset_tmp);
                --index_x;
            }
            for (DefLUint iy = 0; iy < i_max2; ++iy) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindYNeg(bitset_tmp);
                --index_y;
            }
        }
        // top surface
        bitset_tmp2 =  FindZPos(bitset_tmp2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        ++index_z;
        for (DefLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max2; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
        for (DefLUint iy = 0;  iy < i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                ++index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            index_x = indices[kXIndex];
            ++index_y;
        }
    }
}
/**
 * @brief function to check if a node is on the interface of partitioned blocks at the current rank
 * @param[in] sfbitset_in space filling code of input node.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @return true if the input node is on the interface of partitioned blocks at the current rank
 */
bool SFBitsetAux3D::CheckPartitionLimits(const DefSFBitset& sfbitset_in,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefLUint, 3>& code_domain_min, const std::array<DefLUint, 3>& code_domain_max) const {
    std::array<DefSFBitset, 27> neighbours;
    DefSFCodeToUint code_neighbour;
    std::array<DefLUint, 3> indices;
    SFBitsetFindAllNeighbours(sfbitset_in, &neighbours);

    for (DefUint i = 1; i < 27; ++i) {
        code_neighbour = neighbours[i].to_ullong();
        if (code_neighbour < code_partition_min || code_neighbour > code_partition_max) {
            SFBitsetComputeIndices(neighbours[i], &indices);
            if (indices[kXIndex] < code_domain_min[kXIndex] || indices[kXIndex] > code_domain_max[kXIndex]
                || indices[kYIndex] < code_domain_min[kYIndex] || indices[kYIndex] > code_domain_max[kYIndex]
                || indices[kZIndex] < code_domain_min[kZIndex] || indices[kZIndex] > code_domain_max[kZIndex]) {
            } else {
                return true;
            }
        }
    }
    return false;
}

#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
