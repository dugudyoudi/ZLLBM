//  Copyright (c) 2022, Zhengliang Liu
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
#include "auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
#ifdef DEBUG_CHECK_GRID
#include "io/log_write.h"
#endif  // DEBUG_CHECK_GRID
namespace rootproject {
namespace amrproject {
namespace grid {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief initialize morton code related varialbes.
*/
void SFBitsetAux2D::SFBitsetInitial() {
    // reference morton code used to take digitals at a given direction
    k0SFBitsetTakeXRef_.at(kRefCurrent_) = 0;
    k0SFBitsetTakeYRef_.at(kRefCurrent_) = 0;

    for (DefSizet i = 0; i < kSFBitsetBit / 2; ++i) {
        DefSizet pos = i * 2;
        k0SFBitsetTakeXRef_.at(kRefCurrent_).set(pos, true);
        k0SFBitsetTakeYRef_.at(kRefCurrent_).set(pos + 1, true);
    }
    k0SfBItsetTakeCenter_ = DefSFBitset(std::string("0011"));
    k0SFBitsetTakeXRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeXRef_.at(kRefCurrent_);
    k0SFBitsetTakeYRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeYRef_.at(kRefCurrent_);
}
/**
* @brief function to generate morton code based 2D indices.
* @param[in]  coordi_index      indices related to cooridinates.
* @return  morton code.
*/
DefSFBitset SFBitsetAux2D::SFBitsetEncoding(
    const std::array<DefLUint, 2>& coordi_index) {
    DefSFBitset sfbitset_code = 0;
    for (DefSizet i = 0; i < (kSFBitsetBit / 2); ++i) {
        sfbitset_code |= ((coordi_index.at(kXIndex) &
            ((static_cast<DefLUint>(1)) << i)) << i)
            | ((coordi_index.at(kYIndex) &
                ((static_cast<DefLUint>(1)) << i)) << (i + 1));
    }
    return sfbitset_code;
}
/**
* @brief function to compute the 2D coordinates from morton code.
* @param[in]  bitset_in           morton code of current node.
* @param[in]  grid_space    grid spacing in x, y and z direction
*                           at a given refinement level.
* @param[out]  ptr_vec_coordi     pointer to 2D coordinates.
*/
void SFBitsetAux2D::SFBitsetComputeCoordinate(const DefSFBitset& bitset_in,
    const std::array<DefReal, 2>& grid_space,
    std::array<DefReal, 2>* const ptr_coordi) {
    *ptr_coordi = { 0., 0. };
    for (DefSizet i = 0; i < kSFBitsetBit / 2; ++i) {
        if (bitset_in.test(i * 2)) {
            ptr_coordi->at(kXIndex) += grid_space.at(kXIndex)
                * TwoPowerN(i);
        }
        if (bitset_in.test(i * 2 + 1)) {
            ptr_coordi->at(kYIndex) += grid_space.at(kYIndex)
                * TwoPowerN(i);
        }
    }
}
/**
* @brief function to compute the 2D coordinates from morton code.
* @param[in]  bitset_in         morton code of current node.
* @param[out]  ptr_vec_indices   pointer to 2D Indices.
*/
void SFBitsetAux2D::SFBitsetComputeIndices(
    const DefSFBitset& bitset_in,
    std::array<DefLUint, 2>* const ptr_indices) {
    *ptr_indices = { 0, 0 };
    for (DefSizet i = 0; i < kSFBitsetBit / 2; ++i) {
        if (bitset_in.test(i * 2)) {
            ptr_indices->at(kXIndex) +=
                TwoPowerN(static_cast<DefLUint>(i));
        }
        if (bitset_in.test(i * 2 + 1)) {
            ptr_indices->at(kYIndex) +=
                TwoPowerN(static_cast<DefLUint>(i));
        }
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
void SFBitsetAux3D::SFBitsetInitial() {
    k0SFBitsetTakeXRef_.at(kRefCurrent_) = 0;
    k0SFBitsetTakeYRef_.at(kRefCurrent_) = 0;
    k0SFBitsetTakeZRef_.at(kRefCurrent_) = 0;
    for (DefSizet i = 0; i < kSFBitsetBit / 3; ++i) {
        DefSizet pos = i * 3;
        k0SFBitsetTakeXRef_.at(kRefCurrent_).set(pos, true);
        k0SFBitsetTakeYRef_.at(kRefCurrent_).set(pos + 1, true);
        k0SFBitsetTakeZRef_.at(kRefCurrent_).set(pos + 2, true);
    }
    k0SfBItsetTakeCenter_ = DefSFBitset(std::string("0111"));

    k0SFBitsetTakeXRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeXRef_.at(kRefCurrent_);
    k0SFBitsetTakeYRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeYRef_.at(kRefCurrent_);
    k0SFBitsetTakeZRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeZRef_.at(kRefCurrent_);
}
/**
* @brief function to generate morton code based 3D indices.
* @param[in]  coordi_index      indices related to cooridinates.
* @return  morton code.
*/
DefSFBitset SFBitsetAux3D::SFBitsetEncoding(
    const std::array<DefLUint, 3>& coordi_index) {
    DefSFBitset sfbitset_code = 0;
    for (DefSizet i = 0; i < (kSFBitsetBit / 3); ++i) {
        sfbitset_code |= ((coordi_index.at(kXIndex) &
            ((static_cast<DefLUint>(1)) << i)) << (2 * i)) |
            ((coordi_index.at(kYIndex) &
                ((static_cast<DefLUint>(1)) << i)) << (2 * i + 1)) |
            ((coordi_index.at(kZIndex) &
                ((static_cast<DefLUint>(1)) << i)) << (2 * i + 2));
    }
    return sfbitset_code;
}
/**
* @brief function to compute the 3D coordinates from morton code.
* @param[in]  bitset_in          morton code of current node.
* @param[in]  grid_space   grid spacing in x, y and z direction
*                          at a given refinement level.
* @param[out] ptr_vec_coordi  pointer to 3D coordinates.
*/
void SFBitsetAux3D::SFBitsetComputeCoordinate(const DefSFBitset& bitset_in,
     const std::array<DefReal, 3>& grid_space,
    std::array<DefReal, 3>* const ptr_coordi) {
    *ptr_coordi = { 0., 0., 0.};
    for (DefSizet i = 0; i < kSFBitsetBit / 3; ++i) {
        if (bitset_in.test(i * 3)) {
            ptr_coordi->at(kXIndex) +=
                grid_space.at(kXIndex) * TwoPowerN(i);
        }
        if (bitset_in.test(i * 3 + 1)) {
            ptr_coordi->at(kYIndex) +=
                grid_space.at(kYIndex) * TwoPowerN(i);
        }
        if (bitset_in.test(i * 3 + 2)) {
            ptr_coordi->at(kZIndex) +=
                grid_space.at(kZIndex) * TwoPowerN(i);
        }
    }
}
/**
* @brief function to compute the 3D coordinates from morton code.
* @param[in]  bitset_in        morton code of current node.
* @param[out] Indices    pointer to 3D Indices.
*/
void SFBitsetAux3D::SFBitsetComputeIndices(
    const DefSFBitset& bitset_in,
    std::array<DefLUint, 3>* const ptr_indices) {
    *ptr_indices = { 0, 0, 0 };
    for (DefSizet i = 0; i < kSFBitsetBit / 3; ++i) {
        if (bitset_in.test(i * 3)) {
            ptr_indices->at(kXIndex) +=
                TwoPowerN(static_cast<DefLUint>(i));
        }
        if (bitset_in.test(i * 3 + 1)) {
            ptr_indices->at(kYIndex) +=
                TwoPowerN(static_cast<DefLUint>(i));
        }
        if (bitset_in.test(i * 3 + 2)) {
            ptr_indices->at(kZIndex) +=
                TwoPowerN(static_cast<DefLUint>(i));
        }
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namsapce grid
}  // end namespace amrproject
}  // end namespace rootproject
