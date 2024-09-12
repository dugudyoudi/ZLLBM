//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.cpp
* @author Zhengliang Liu
* @brief functions used to manipulate morton codes.
* @date  2022-5-16
*/
#include <string>
#include <array>
#include "auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief initialize morton code related variables.
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
DefSFBitset SFBitsetAux2D::SFBitsetBitsForRefinement(const DefInt i_level) const {
    DefSFBitset bitset = 0;
    for (DefInt i = 0; i < 2 * i_level; ++i) {
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
    const std::array<DefAmrLUint, 2>& coordi_index) const {
    DefSFBitset sfbitset_code = 0;
    for (auto i = 0; i < (kSFBitsetBit / 2); ++i) {
        sfbitset_code |= ((coordi_index.at(kXIndex)
            &((static_cast<DefAmrLUint>(1)) << i)) << i)
            |((coordi_index.at(kYIndex)
            &((static_cast<DefAmrLUint>(1)) << i)) << (i + 1));
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
    std::array<DefAmrLUint, 2> coordi_index =
        {static_cast<DefAmrLUint>(coordi.at(kXIndex) / grid_space[kXIndex] + kEps),
        static_cast<DefAmrLUint>(coordi.at(kYIndex) / grid_space[kYIndex] + kEps)};
    return this->SFBitsetEncoding(coordi_index);
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
    std::array<DefAmrLUint, 2>* const ptr_indices) const {
    *ptr_indices = { 0, 0 };
    for (DefSizet i = 0; i < kSFBitsetBit / 2; ++i) {
        if (bitset_in.test(i * 2)) {
            ptr_indices->at(kXIndex) += static_cast<DefAmrLUint>(TwoPowerN(i));
        }
        if (bitset_in.test(i * 2 + 1)) {
            ptr_indices->at(kYIndex) += static_cast<DefAmrLUint>(TwoPowerN(i));
        }
    }
}
/**
 * @brief function to reset the indices exceeding the computational domain.
 * @param[in] domain_min_indices minimum indices of the computational domain.
 * @param[in] domain_max_indices maximum indices of the computational domain.   
 * @param[out] ptr_i_code pointer to a number to record increment of the space filling code
 * @param[out] ptr_bitset_tmp pointer to space filling code
 * @return an integer indicating the success of the operation
 */
int SFBitsetAux2D::ResetIndicesExceedingDomain(const std::array<DefAmrLUint, 2>& domain_min_indices,
    const std::array<DefAmrLUint, 2>& domain_max_indices,
    DefSFCodeToUint* const ptr_i_code, DefSFBitset* ptr_bitset_tmp) const {
    std::array<DefAmrLUint, 2> indices;
    DefSFCodeToUint& i_code = *ptr_i_code;
    DefSFBitset& sfbitset_tmp = *ptr_bitset_tmp;
    sfbitset_tmp = static_cast<DefSFBitset>(i_code);
    DefInt iter_count;
    DefInt sfcode_block;  // block size of space filling code
    bool bool_exceed_x, bool_exceed_y;
    SFBitsetComputeIndices(sfbitset_tmp, &indices);
    bool_exceed_x =
        indices[kXIndex] > domain_max_indices[kXIndex];
    bool_exceed_y =
        indices[kYIndex] > domain_max_indices[kYIndex];
    iter_count = 0;
    while ((bool_exceed_x || bool_exceed_y)
        && iter_count < max_reset_code_) {
        if (bool_exceed_x) {
            sfcode_block = 4;
            while ((i_code % sfcode_block) == 0) {
                sfcode_block *= 4;
            }
            i_code += sfcode_block / 4;
            sfbitset_tmp = static_cast<DefSFBitset>(i_code);
            SFBitsetComputeIndices(sfbitset_tmp, &indices);
            bool_exceed_x =
                indices[kXIndex] > domain_max_indices[kXIndex];
            bool_exceed_y =
                indices[kYIndex] > domain_max_indices[kYIndex];
        }
        if (bool_exceed_y) {
            sfcode_block = 4;
            while ((i_code % sfcode_block) == 0) {
                sfcode_block *= 4;
            }
            i_code += sfcode_block / 4 * 2;
            sfbitset_tmp = static_cast<DefSFBitset>(i_code);
            SFBitsetComputeIndices(sfbitset_tmp, &indices);
            bool_exceed_x =
                indices[kXIndex] > domain_max_indices[kXIndex];
            bool_exceed_y =
                indices[kYIndex] > domain_max_indices[kYIndex];
        }
        ++iter_count;
    }
#ifdef DEBUG_CHECK_GRID
    if (iter_count > max_reset_code_) {
        return -1;
    }
#endif
    if (indices[kXIndex] < domain_min_indices[kXIndex]) {
        sfbitset_tmp = SFBitsetEncoding(
            { domain_min_indices[kXIndex], indices[kYIndex] });
        i_code = sfbitset_tmp.to_ullong();
        SFBitsetComputeIndices(sfbitset_tmp, &indices);
    }
    if (indices[kYIndex] < domain_min_indices[kYIndex]) {
        sfbitset_tmp = SFBitsetEncoding(
            { indices[kXIndex] , domain_min_indices[kYIndex] });
        i_code = sfbitset_tmp.to_ullong();
        SFBitsetComputeIndices(sfbitset_tmp, &indices);
    }
    return 0;
}
/**
 * @brief function to find coarse nodes lined to the fine nodes.
 * @param[in] sfbitset_fine space filling code of the fine nodes.
 * @param[out] ptr_vec_coarse pointer to space filling codes of the coarse nodes.
 */
void SFBitsetAux2D::FindNeighboringCoarseFromFine(const DefSFBitset& sfbitset_fine,
    std::vector<DefSFBitset>* const ptr_vec_coarse) const {
    ptr_vec_coarse->clear();
    bool bool_x = sfbitset_fine.test(0);
    bool bool_y = sfbitset_fine.test(1);
    if (bool_x || bool_y) {
        DefSFBitset sfbitset_coarse = SFBitsetToOneLowerLevel(sfbitset_fine);
        ptr_vec_coarse->push_back(sfbitset_coarse);
        if (bool_x) {
            DefSFBitset sfbitset_coarse_tmp = FindXPos(sfbitset_coarse);
            ptr_vec_coarse->push_back(sfbitset_coarse_tmp);
            if (bool_y) {
                ptr_vec_coarse->push_back(FindYPos(sfbitset_coarse_tmp));
            }
        }
        if (bool_y) {
            ptr_vec_coarse->push_back(FindYPos(sfbitset_coarse));
        }
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief initialize morton code related variables.
*/
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
DefSFBitset SFBitsetAux3D::SFBitsetBitsForRefinement(const DefInt max_level) const {
    DefSFBitset bitset = 0;
    for (DefInt i = 0; i < 3 * max_level; ++i) {
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
    const std::array<DefAmrLUint, 3>& coordi_index)  const {
    DefSFBitset sfbitset_code = 0;
    for (DefSizet i = 0; i < (kSFBitsetBit / 3); ++i) {
        sfbitset_code |= ((coordi_index.at(kXIndex) &
            ((static_cast<DefAmrLUint>(1)) << i)) << (2 * i)) |
            ((coordi_index.at(kYIndex) &
                ((static_cast<DefAmrLUint>(1)) << i)) << (2 * i + 1)) |
            ((coordi_index.at(kZIndex) &
                ((static_cast<DefAmrLUint>(1)) << i)) << (2 * i + 2));
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
    std::array<DefAmrLUint, 3> coordi_index =
    { static_cast<DefAmrLUint>(coordi.at(kXIndex) / grid_space[kXIndex] + kEps),
      static_cast<DefAmrLUint>(coordi.at(kYIndex) / grid_space[kYIndex] + kEps),
      static_cast<DefAmrLUint>(coordi.at(kZIndex) / grid_space[kZIndex] + kEps) };
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
    std::array<DefAmrLUint, 3>* const ptr_indices) const {
    *ptr_indices = { 0, 0, 0 };
    for (DefSizet i = 0; i < kSFBitsetBit / 3; ++i) {
        if (bitset_in.test(i * 3)) {
            ptr_indices->at(kXIndex) += static_cast<DefAmrLUint>(TwoPowerN(i));
        }
        if (bitset_in.test(i * 3 + 1)) {
            ptr_indices->at(kYIndex) += static_cast<DefAmrLUint>(TwoPowerN(i));
        }
        if (bitset_in.test(i * 3 + 2)) {
            ptr_indices->at(kZIndex) += static_cast<DefAmrLUint>(TwoPowerN(i));
        }
    }
}
/**
 * @brief function to reset the indices exceeding the computational domain.
* @param[in] domain_min_indices minimum indices of the computational domain.
 * @param[in] domain_max_indices maximum indices of the computational domain.  
 * @param ptr_i_code pointer to a number to record increment of the space filling code
 * @param ptr_bitset_tmp pointer to space filling code
 * @return an integer indicating the success of the operation
 */
int SFBitsetAux3D::ResetIndicesExceedingDomain(const std::array<DefAmrLUint, 3>& domain_min_indices,
    const std::array<DefAmrLUint, 3>& domain_max_indices,
    DefSFCodeToUint* const ptr_i_code, DefSFBitset* ptr_bitset_tmp) const {
    std::array<DefAmrLUint, 3> indices;
    DefSFCodeToUint& i_code = *ptr_i_code;
    DefSFBitset& sfbitset_tmp = *ptr_bitset_tmp;
    sfbitset_tmp = static_cast<DefSFBitset>(i_code);
    DefInt iter_count;
    DefInt sfcode_block;  // block size of space filling code
    bool bool_exceed_x, bool_exceed_y, bool_exceed_z;
    SFBitsetComputeIndices(sfbitset_tmp, &indices);
    bool_exceed_x =
        indices[kXIndex] > domain_max_indices[kXIndex];
    bool_exceed_y =
        indices[kYIndex] > domain_max_indices[kYIndex];
    bool_exceed_z =
        indices[kZIndex] > domain_max_indices[kZIndex];
    iter_count = 0;
    while ((bool_exceed_x || bool_exceed_y || bool_exceed_z)
        && iter_count < max_reset_code_) {
        if (bool_exceed_x) {
            sfcode_block = 8;
            while ((i_code % sfcode_block) == 0) {
                sfcode_block *= 8;
            }
            i_code += sfcode_block / 8;
            sfbitset_tmp = static_cast<DefSFBitset>(i_code);
            SFBitsetComputeIndices(sfbitset_tmp, &indices);
            bool_exceed_x =
                indices[kXIndex] > domain_max_indices[kXIndex];
            bool_exceed_y =
                indices[kYIndex] > domain_max_indices[kYIndex];
            bool_exceed_z =
                indices[kZIndex] > domain_max_indices[kZIndex];
        }
        if (bool_exceed_y) {
            sfcode_block = 8;
            while ((i_code % sfcode_block) == 0) {
                sfcode_block *= 8;
            }
            i_code += sfcode_block / 8 * 2;
            sfbitset_tmp = static_cast<DefSFBitset>(i_code);
            SFBitsetComputeIndices(sfbitset_tmp, &indices);
            bool_exceed_x =
                indices[kXIndex] > domain_max_indices[kXIndex];
            bool_exceed_y =
                indices[kYIndex] > domain_max_indices[kYIndex];
            bool_exceed_z =
                indices[kZIndex] > domain_max_indices[kZIndex];
        }
        if (bool_exceed_z) {
            sfcode_block = 8;
            while ((i_code % sfcode_block) == 0) {
                sfcode_block *= 8;
            }
            i_code += sfcode_block / 8 * 4;
            sfbitset_tmp = static_cast<DefSFBitset>(i_code);
            SFBitsetComputeIndices(sfbitset_tmp, &indices);
            bool_exceed_x =
                indices[kXIndex] > domain_max_indices[kXIndex];
            bool_exceed_y =
                indices[kYIndex] > domain_max_indices[kYIndex];
            bool_exceed_z =
                indices[kZIndex] > domain_max_indices[kZIndex];
        }
        ++iter_count;
    }
#ifdef DEBUG_CHECK_GRID
    if (iter_count > max_reset_code_) {
        return -1;
    }
#endif
    if (indices[kXIndex] < domain_min_indices[kXIndex]) {
        sfbitset_tmp = SFBitsetEncoding(
            { domain_min_indices[kXIndex], indices[kYIndex], indices[kZIndex] });
        i_code = sfbitset_tmp.to_ullong();
        SFBitsetComputeIndices(sfbitset_tmp, &indices);
    }
    if (indices[kYIndex] < domain_min_indices[kYIndex]) {
        sfbitset_tmp = SFBitsetEncoding(
            { indices[kXIndex], domain_min_indices[kYIndex], indices[kZIndex] });
        i_code = sfbitset_tmp.to_ullong();
        SFBitsetComputeIndices(sfbitset_tmp, &indices);
    }
    if (indices[kZIndex] < domain_min_indices[kZIndex]) {
        sfbitset_tmp = SFBitsetEncoding(
            { indices[kXIndex], indices[kYIndex], domain_min_indices[kZIndex] });
        i_code = sfbitset_tmp.to_ullong();
        SFBitsetComputeIndices(sfbitset_tmp, &indices);
    }
    return 0;
}
/**
 * @brief function to find coarse nodes lined to the fine nodes.
 * @param[in] sfbitset_fine space filling code of the fine nodes.
 * @param[out] ptr_vec_coarse pointer to space filling codes of the coarse nodes.
 */
void SFBitsetAux3D::FindNeighboringCoarseFromFine(const DefSFBitset& sfbitset_fine,
    std::vector<DefSFBitset>* const ptr_vec_coarse) const {
    ptr_vec_coarse->clear();
    bool bool_x = sfbitset_fine.test(0);
    bool bool_y = sfbitset_fine.test(1);
    bool bool_z = sfbitset_fine.test(2);
    if (bool_x || bool_y || bool_z) {
        DefSFBitset sfbitset_coarse = SFBitsetToOneLowerLevel(sfbitset_fine);
        ptr_vec_coarse->push_back(sfbitset_coarse);
        if (bool_x) {
            DefSFBitset sfbitset_coarse_tmp = FindXPos(sfbitset_coarse);
            ptr_vec_coarse->push_back(sfbitset_coarse_tmp);
            if (bool_y) {
                DefSFBitset sfbitset_coarse_tmp1 = FindYPos(sfbitset_coarse_tmp);
                ptr_vec_coarse->push_back(sfbitset_coarse_tmp1);
                if (bool_z) {
                    ptr_vec_coarse->push_back(FindZPos(sfbitset_coarse_tmp1));
                }
            }
        }
        if (bool_y) {
            DefSFBitset sfbitset_coarse_tmp = FindYPos(sfbitset_coarse);
            ptr_vec_coarse->push_back(FindYPos(sfbitset_coarse));
            if (bool_z) {
                ptr_vec_coarse->push_back(FindZPos(sfbitset_coarse_tmp));
            }
        }
        if (bool_z) {
            ptr_vec_coarse->push_back(FindZPos(sfbitset_coarse));
        }
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
