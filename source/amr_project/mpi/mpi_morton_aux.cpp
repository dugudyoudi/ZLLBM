//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_morton_aux.cpp
* @author Zhengliang Liu
* @brief functions used to manipulate morton codes for MPI partition.
* @date  2023-7-16
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
void SFBitsetAux2D::FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefAmrIndexLUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length - 1;
    code_remain = (code_in - 1) % 4;
    std::array<DefAmrIndexLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 0:
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1 + 1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1 + 1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
void SFBitsetAux2D::FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefAmrIndexLUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length;
    code_remain = code_in % 4;
    std::array<DefAmrIndexLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    bitset_tmp = bitset_current;
    switch (code_remain) {
        case 3:
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1 + 1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1 + 1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
void SFBitsetAux2D::FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefAmrIndexLUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
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
void SFBitsetAux2D::FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefAmrIndexLUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_current, bitset_tmp;
    DefSFCodeToUint code_current;
    code_current = code_in * block_length * block_length;
    code_remain = code_in % 4;
    std::array<DefAmrIndexLUint, 2> indices;
    bitset_current = static_cast<DefSFBitset>(code_current);
    SFBitsetComputeIndices(bitset_current, &indices);
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
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
 * @param[in] bitset_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_length length of the current block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux2D::FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefAmrIndexLUint block_length,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
         std::array<DefAmrIndexLUint, 2> indices;
    DefSFBitset bitset_tmp = bitset_corner;
    SFBitsetComputeIndices(bitset_corner, &indices);
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
        && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
        && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
            code_domain_min, code_domain_max)
        && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
        ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        DefAmrIndexLUint i_max1 = block_length - 1;
        for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
        for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
        for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
        for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
    const std::array<DefAmrIndexLUint, 2>& code_domain_min, const std::array<DefAmrIndexLUint, 2>& code_domain_max) const {
    std::array<DefSFBitset, 9> neighbors;
    DefSFCodeToUint code_neighbour;
    std::array<DefAmrIndexLUint, 2> indices;
    SFBitsetFindAllNeighbors(sfbitset_in, &neighbors);
    for (DefAmrIndexUint i = 1; i < 9; ++i) {
        code_neighbour = neighbors[i].to_ullong();
        if (code_neighbour < code_partition_min || code_neighbour > code_partition_max) {
            SFBitsetComputeIndices(neighbors[i], &indices);
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
/**
 * @brief function to find interface of partitioned blocks by several blocks at given level based on the maximum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFBitset bitset_tmp;
    DefSFCodeToUint code_current;
    const DefAmrIndexLUint block_length = 1 << block_level;
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
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    const DefAmrIndexLUint block_length = 1 << block_level;
    const DefAmrIndexLUint block_length_cubic =  block_length * block_length * block_length;
    code_remain = code_in % 8;
    code_current = (code_in + 8 - code_remain)* block_length_cubic - 1;
    DefSFBitset bitset_tmp = static_cast<DefSFBitset>(code_current);

    std::array<DefAmrIndexLUint, 3> indices;
    SFBitsetComputeIndices(bitset_tmp, &indices);
    switch (code_remain) {
        case 7:
            bitset_tmp = static_cast<DefSFBitset>(code_in * block_length_cubic);
            FindPartitionOneBlock(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 6:
            FindPartition2BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 5:
            FindPartition3BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 4:
            FindPartition4BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 3:
            FindPartition5BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 2:
            FindPartition6BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 1:
            FindPartition7BlocksMin(bitset_tmp, block_level, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        case 0:
            bitset_tmp = static_cast<DefSFBitset>(code_in * block_length_cubic);
            FindPartitionOneBlock(bitset_tmp, block_level + 1, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max, ptr_map_partitioned);
            break;
        default:
            break;
    }
}
/**
 * @brief function to find interface of partitioned blocks one by one block at given level based on the maximum space fill code. 
 * @param[in] code_in space filling code of the input node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    const DefAmrIndexLUint block_length = 1 << block_level;
    const DefAmrIndexLUint block_length_cubic =  block_length * block_length * block_length;
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
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    DefSFCodeToUint code_remain;
    DefSFCodeToUint code_current;
    const DefAmrIndexLUint block_length = 1 << block_level;
    const DefAmrIndexLUint block_length_cubic =  block_length * block_length * block_length;
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
 * @param[in] bitset_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1;
    SFBitsetComputeIndices(bitset_corner, &indices);
    DefSFBitset bitset_tmp = bitset_corner, bitset_tmp2 = bitset_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
 * @brief function to find interface of partitioned blocks by two blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition2BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
 * @brief function to find interface of partitioned blocks by three blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition3BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        index_x = indices[kXIndex];
        ++index_y;
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix <= i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy <= i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max2; ++iy) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
 * @brief function to find interface of partitioned blocks by four blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition4BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        index_x = indices[kXIndex];
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
 * @brief function to find interface of partitioned blocks by five blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition5BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        ++index_z;
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            ++index_y;
            index_x = indices[kXIndex];
        }
        // lower vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_min_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            ++index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindXPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            ++index_y;
            index_x = indices2[kXIndex];
        }
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            ++index_y;
            index_x = indices2[kXIndex];
        }
        // upper vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            bitset_tmp2 = FindZPos(bitset_tmp2);
            ++index_z;
        }
        // top surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            ++index_y;
            index_x = indices[kXIndex];
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by six blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition6BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        ++index_z;
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            ++index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            ++index_y;
            index_x = indices2[kXIndex];
        }
         // upper vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            ++index_z;
        }
        // top surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            ++index_y;
            index_x = indices[kXIndex];
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by seven blocks with increasing space filling code. 
 * @param[in] bitset_min_corner space filling code of lower left (0, 0, 0) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition7BlocksMax(const DefSFBitset& bitset_min_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    SFBitsetComputeIndices(bitset_min_corner, &indices);
    DefSFBitset bitset_tmp = bitset_min_corner, bitset_tmp2 = bitset_min_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
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
        ++index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYPos(bitset_min_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZPos(bitset_min_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        ++index_z;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp2 = FindXPos(bitset_tmp);
        ++index_x;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp2, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp2) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp2, 0});
        }
        bitset_tmp2 = FindYPos(bitset_tmp);
        ++index_y;
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            ++index_z;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_min_corner);
        bitset_tmp = FindXPos(bitset_tmp);
        bitset_tmp = FindYPos(bitset_tmp);
        bitset_tmp = FindZPos(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
        bitset_tmp2 = FindZNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix <= i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy <= i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max2; ++iy) {
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
            bitset_tmp2 =  FindZPos(bitset_tmp2);
            ++index_z;
        }
        // top surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
 * @brief function to find interface of partitioned blocks by two blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition2BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
        }
        // bottom surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZNeg(bitset_tmp2);
        --index_z;
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by three blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition3BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        --index_y;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            bitset_tmp = bitset_tmp2;
            --index_z;
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix <= i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy <= i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max2; ++iy) {
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
        }
        // bottom surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZNeg(bitset_tmp2);
        --index_z;
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by four blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition4BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        // vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        }
        // bottom surface
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 =  FindZNeg(bitset_tmp2);
        --index_z;
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by five blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition5BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
         if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        --index_z;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            --index_y;
            index_x = indices[kXIndex];
        }
        // upper vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        }
        // middle bottom surface
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
                if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
                    && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
                    && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
                    && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                        code_domain_min, code_domain_max)
                    && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
                    ptr_map_partitioned->insert({bitset_tmp, 0});
                }
                bitset_tmp = FindXPos(bitset_tmp);
                --index_x;
            }
            bitset_tmp2 = FindYPos(bitset_tmp2);
            ++index_y;
            index_x = indices2[kXIndex];
        }
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindYNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            --index_y;
            index_x = indices2[kXIndex];
        }
        // lower vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            bitset_tmp2 = FindZNeg(bitset_tmp2);
            --index_z;
        }
        // bottom surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            ++index_y;
            index_x = indices[kXIndex];
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by six blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition6BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length *2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        --index_z;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            --index_y;
            index_x = indices[kXIndex];
        }
        // upper vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        }
        // middle bottom surface
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindYNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            --index_y;
            index_x = indices2[kXIndex];
        }
         // lower vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max1; ++iy) {
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
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
        }
        // bottom surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            ++index_y;
            index_x = indices[kXIndex];
        }
    }
}
/**
 * @brief function to find interface of partitioned blocks by seven blocks with decreasing space filling code. 
 * @param[in] bitset_max_corner space filling code of upper right (x, y, z) conner node.
 * @param[in] block_level level of the block.
 * @param[in] code_partition_min minimum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_partition_max maximum space fill code of the partitioned blocks of the current rank.
 * @param[in] code_domain_min minimum indicies of computational domain.
 * @param[in] code_domain_max maximum indicies of computational domain.
 * @param[out] ptr_map_partitioned  map of unsigned integers to unsigned integers representing the partitioned results.
 */
void SFBitsetAux3D::FindPartition7BlocksMin(const DefSFBitset& bitset_max_corner, const DefAmrIndexUint block_level,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max,
    DefMap<DefAmrIndexUint>* const ptr_map_partitioned) const {
    std::array<DefAmrIndexLUint, 3> indices;
    const DefAmrIndexLUint block_length = 1 << block_level;
    DefAmrIndexLUint i_max1 = block_length - 1, i_max2 = block_length * 2 - 1;
    SFBitsetComputeIndices(bitset_max_corner, &indices);
    DefSFBitset bitset_tmp = bitset_max_corner, bitset_tmp2 = bitset_max_corner;
    DefAmrIndexLUint index_x = indices[kXIndex], index_y = indices[kYIndex], index_z = indices[kZIndex];
    if (block_length == 1) {
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindXNeg(bitset_max_corner);
        --index_x;
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
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindYNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp = FindZNeg(bitset_max_corner);
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        --index_z;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp, 0});
        }
        bitset_tmp2 = FindXNeg(bitset_tmp);
        --index_x;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp2, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp2) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp2, 0});
        }
        bitset_tmp2 = FindYNeg(bitset_tmp);
        --index_y;
        if (index_x >= code_domain_min[kXIndex] && index_x <= code_domain_max[kXIndex]
            && index_y >= code_domain_min[kYIndex] && index_y <= code_domain_max[kYIndex]
            && index_z >= code_domain_min[kZIndex] && index_z <= code_domain_max[kZIndex]
            && CheckPartitionLimits(bitset_tmp2, code_partition_min, code_partition_max,
                code_domain_min, code_domain_max)
            && ptr_map_partitioned->find(bitset_tmp2) == ptr_map_partitioned->end()) {
            ptr_map_partitioned->insert({bitset_tmp2, 0});
        }
    } else {
        // top surface
        for (DefAmrIndexLUint iy = 0;  iy <= i_max2; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices[kXIndex];
            --index_y;
        }
        // upper vertical surfaces
        index_x = indices[kXIndex];
        index_y = indices[kYIndex];
        bitset_tmp2 = bitset_max_corner;
        for (DefAmrIndexLUint iz = 1; iz < i_max1; ++iz) {
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0;  ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0;  iy < i_max2; ++iy) {
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
        }
        // middle bottom surface
        std::array<DefAmrIndexLUint, 3> indices2;
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindYNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices2[kXIndex];
            --index_y;
        }
        // lower vertical surfaces
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindXNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        index_z = indices2[kZIndex];
        for (DefAmrIndexLUint iz = 0; iz <= i_max1; ++iz) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint iy = 0; iy < i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max2; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy < i_max2; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix < i_max1; ++ix) {
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
            for (DefAmrIndexLUint iy = 0; iy <= i_max1; ++iy) {
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
            for (DefAmrIndexLUint ix = 0; ix <= i_max1; ++ix) {
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
            bitset_tmp2 =  FindZNeg(bitset_tmp2);
            --index_z;
        }
        // bottom surface
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max2; ++ix) {
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
        bitset_tmp = SFBitsetToNLowerLevel(block_level, bitset_max_corner);
        bitset_tmp = FindZNeg(bitset_tmp);
        bitset_tmp2 = SFBitsetToNHigherLevel(block_level, bitset_tmp);
        bitset_tmp2 = FindYNeg(bitset_tmp2);
        SFBitsetComputeIndices(bitset_tmp2, &indices2);
        index_x = indices2[kXIndex];
        index_y = indices2[kYIndex];
        for (DefAmrIndexLUint iy = 0;  iy <= i_max1; ++iy) {
            bitset_tmp = bitset_tmp2;
            for (DefAmrIndexLUint ix = 0;  ix <= i_max1; ++ix) {
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
            bitset_tmp2 = FindYNeg(bitset_tmp2);
            index_x = indices2[kXIndex];
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
bool SFBitsetAux3D::CheckPartitionLimits(const DefSFBitset& sfbitset_in,
    const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
    const std::array<DefAmrIndexLUint, 3>& code_domain_min, const std::array<DefAmrIndexLUint, 3>& code_domain_max) const {
    std::array<DefSFBitset, 27> neighbors;
    DefSFCodeToUint code_neighbour;
    std::array<DefAmrIndexLUint, 3> indices;
    SFBitsetFindAllNeighbors(sfbitset_in, &neighbors);

    for (DefAmrIndexUint i = 1; i < 27; ++i) {
        code_neighbour = neighbors[i].to_ullong();
        if (code_neighbour < code_partition_min || code_neighbour > code_partition_max) {
            SFBitsetComputeIndices(neighbors[i], &indices);
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
