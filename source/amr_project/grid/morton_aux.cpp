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
namespace rootproject {
namespace amrproject {
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
    k0SfBitsetCurrentLevelBits_ = DefSFBitset(std::string("0011"));
    k0SFBitsetTakeXRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeXRef_.at(kRefCurrent_);
    k0SFBitsetTakeYRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeYRef_.at(kRefCurrent_);
}
/**
* @brief take bits representing coordinates at refinement levels 
*        higher level 0(the background level).
* @param[in] i_level level of refinement
*/
DefSFBitset SFBitsetAux2D::SFBitsetBitsForRefinement(
    const DefSizet i_level) const {
    DefSFBitset bitset = 0;
    for (DefSizet i = 0; i < 2 * i_level; ++i) {
        bitset.set(i, true);
    }
    return bitset;
}
/**
* @brief function to generate morton code based 2D indices.
* @param[in]  coordi_index      indices related to coordinates.
* @return  morton code.
*/
DefSFBitset SFBitsetAux2D::SFBitsetEncoding(
    const std::array<DefLUint, 2>& coordi_index) const {
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
* @brief function to generate morton code based 2D coordinates.
* @param[in]  grid_space      grid space of the current grid level.
* @param[in]  coordi      coordinates of a point.
* @param[in]  ptr_bitset_aux     pointer to 2D space filling curve class.
* @return  morton code.
*/
DefSFBitset SFBitsetAux2D::SFBitsetEncodingCoordi(
    const std::vector<DefReal>& grid_space,
    const std::vector<DefReal>& coordi) const {
    std::array<DefLUint, 2> coorid_index =
    { static_cast<DefLUint>(coordi.at(kXIndex) / grid_space[kXIndex] + kEps),
      static_cast<DefLUint>(coordi.at(kYIndex) / grid_space[kYIndex] + kEps)};
    return this->SFBitsetEncoding(coorid_index);
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
    std::array<DefReal, 2>* const ptr_coordi)  const {
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
    std::array<DefLUint, 2>* const ptr_indices) const {
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
    k0SfBitsetCurrentLevelBits_ = DefSFBitset(std::string("0111"));

    k0SFBitsetTakeXRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeXRef_.at(kRefCurrent_);
    k0SFBitsetTakeYRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeYRef_.at(kRefCurrent_);
    k0SFBitsetTakeZRef_.at(kRefOthers_) =
        ~k0SFBitsetTakeZRef_.at(kRefCurrent_);
}
/**
* @brief take bits representing coordinates at refinement levels
*        higher level 0(the background level).
* @param[in] i_level level of refinement
*/
DefSFBitset SFBitsetAux3D::SFBitsetBitsForRefinement(
    const DefSizet max_level) const{
    DefSFBitset bitset = 0;
    for (DefSizet i = 0; i < 3 * max_level; ++i) {
        bitset.set(i, true);
    }
    return bitset;
}
/**
* @brief function to generate morton code based 3D indices.
* @param[in]  coordi_index      indices related to coordinates.
* @return  morton code.
*/
DefSFBitset SFBitsetAux3D::SFBitsetEncoding(
    const std::array<DefLUint, 3>& coordi_index)  const {
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
* @brief function to generate morton code based 3D coordinates.
* @param[in]  grid_space      grid space of the current grid level.
* @param[in]  coordi      coordinates of a point.
* @param[in]  ptr_bitset_aux     pointer to 3D space filling curve class.
* @return  morton code.
*/
DefSFBitset SFBitsetAux3D::SFBitsetEncodingCoordi(
    const std::vector<DefReal>& grid_space,
    const std::vector<DefReal>& coordi) const {
    std::array<DefLUint, 3> coordi_index =
    { static_cast<DefLUint>(coordi.at(kXIndex) / grid_space[kXIndex] + kEps),
      static_cast<DefLUint>(coordi.at(kYIndex) / grid_space[kYIndex] + kEps),
      static_cast<DefLUint>(coordi.at(kZIndex) / grid_space[kZIndex] + kEps) };
    return this->SFBitsetEncoding(coordi_index);
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
    std::array<DefReal, 3>* const ptr_coordi)  const {
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
    std::array<DefLUint, 3>* const ptr_indices) const {
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
}  // end namespace amrproject
}  // end namespace rootproject
