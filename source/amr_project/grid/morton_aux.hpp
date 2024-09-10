//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.hpp
* @author Zhengliang Liu
* @date  2022-8-16
* @brief incline functions to manipulate morton code
*/

#ifndef SOURCE_AMR_PROJECT_GRID_MORTON_AUX_HPP_
#define SOURCE_AMR_PROJECT_GRID_MORTON_AUX_HPP_
#ifdef DEBUG_CHECK_GRID
#include <bit>
#include "io/log_write.h"
#endif  // DEBUG_CHECK_GRID
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to convert morton code to 1 level lower.
* @param[in]  morton_in   input morton code.
* @return   morton code at 1 lower level.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToOneLowerLevel(
    const DefSFBitset& morton_in) const {
    return morton_in >> 2;
}
/**
* @brief function to convert morton code to 1 level higher.
* @param[in]  morton_in   input morton code.
* @return   morton code 1 higher level.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToOneHigherLevel(
    const DefSFBitset& morton_in) const {
    return morton_in << 2;
}
/**.
* @brief function to convert morton code to n level lower.
* @param[in]  n_level  number of levels need to be shrank.
* @param[in]  morton_in   input morton code.
* @return   morton code at n lower levels.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToNLowerLevel(
    const DefInt n_level,
    const DefSFBitset& morton_in)  const {
    return morton_in >> (2 * n_level);
}
/**
* @brief function to convert morton code to n level higher.
* @param[in]  n_level  number of levels need to be stretched.
* @param[in]  morton_in   input morton code.
* @return   morton code at n higher levels.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToNHigherLevel(
    const DefInt n_level, const DefSFBitset& morton_in) const {
    return morton_in << (2 * n_level);
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the negative x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux2D::FindXNeg(
    const DefSFBitset& bitset_in) const {
    DefSFCodeToUint x_ull = (bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_
        .at(SFBitsetAux2D::kRefCurrent_)).to_ullong();
    if (x_ull == 0) {
        return bitset_in;
    } else {
        return ((static_cast<DefSFBitset>
            (x_ull - 1) & SFBitsetAux2D::k0SFBitsetTakeXRef_
            .at(SFBitsetAux2D::kRefCurrent_))
            | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_
                .at(SFBitsetAux2D::kRefCurrent_)));
    }
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the positive x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux2D::FindXPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset sfbitest_tmp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefOthers_)).to_ullong() + 1);
    return ((sfbitest_tmp & SFBitsetAux2D::k0SFBitsetTakeXRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the negative y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux2D::FindYNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFCodeToUint y_ull = (bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_
        .at(SFBitsetAux2D::kRefCurrent_)).to_ullong();
    if (y_ull == 0) {
        return bitset_in;
    } else {
        return ((static_cast<DefSFBitset>
            (y_ull - 1) & SFBitsetAux2D::k0SFBitsetTakeYRef_
            .at(SFBitsetAux2D::kRefCurrent_))
            | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_
                .at(SFBitsetAux2D::kRefCurrent_)));
    }
}
/**
* @brief function to find morton code (2D) of the neighbour in the positive y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux2D::FindYPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset sfbitest_tmp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefOthers_)).to_ullong() + 1);
    return ((sfbitest_tmp & SFBitsetAux2D::k0SFBitsetTakeYRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
/**
* @brief function to find morton code coincide with code at how many lower refinement levels.
* @param[in]  bitset_in   morton code of the current node.
* @return  number of coincident lower refinement levels.
*/
inline DefInt SFBitsetAux2D::SFBitsetCoincideLevel(const DefSFBitset& bitset_in) const {
    return static_cast<DefInt>(std::countr_zero(bitset_in.to_ullong()))/2;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to convert morton code to 1 level lower.
* @param[in]  morton_in   input morton code.
* @return   morton code at 1 lower level.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToOneLowerLevel(
    const DefSFBitset& morton_in) const {
    return morton_in >> 3;
}
/**
* @brief function to convert morton code to 1 level higher.
* @param[in]  morton_in   input morton code.
* @return   morton code at 1 higher level.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToOneHigherLevel(
    const DefSFBitset& morton_in) const {
    return morton_in << 3;
}
/**
* @brief function to convert morton code to n level lower.
* @param[in]  n_level  number of levels need to be shrank.
* @param[in]  morton_in   input morton code.
* @return   morton code at n lower levels.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToNLowerLevel(
    const DefInt n_level, const DefSFBitset& morton_in) const {
    return morton_in >> (3 * n_level);
}
/**
* @brief function to convert morton code to n level higher.
* @param[in]  n_level  number of levels need to be stretched.
* @param[in]  morton_in   input morton code.
* @return   morton code at n higher levels.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToNHigherLevel(
    const DefInt n_level, const DefSFBitset& morton_in)  const {
    return morton_in << (3 * n_level);
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the negative x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindXNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFCodeToUint x_ull = (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_
        .at(SFBitsetAux3D::kRefCurrent_)).to_ullong();
    if (x_ull == 0) {
        return bitset_in;
    } else {
        return ((static_cast<DefSFBitset>
            (x_ull - 1) & SFBitsetAux3D::k0SFBitsetTakeXRef_
            .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_
                .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_
                .at(SFBitsetAux3D::kRefCurrent_)));
    }
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindXPos(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset sfbitest_tmp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((sfbitest_tmp & SFBitsetAux3D::k0SFBitsetTakeXRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the negative y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindYNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFCodeToUint y_ull = (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_
        .at(SFBitsetAux3D::kRefCurrent_)).to_ullong();
    if (y_ull == 0) {
        return bitset_in;
    } else {
        return ((static_cast<DefSFBitset>
            (y_ull - 1) & SFBitsetAux3D::k0SFBitsetTakeYRef_
            .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_
                .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_
                .at(SFBitsetAux3D::kRefCurrent_)));
    }
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindYPos(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset sfbitest_tmp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((sfbitest_tmp & SFBitsetAux3D::k0SFBitsetTakeYRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the negative z direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindZNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFCodeToUint z_ull = (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_
        .at(SFBitsetAux3D::kRefCurrent_)).to_ullong();
    if (z_ull == 0) {
        return bitset_in;
    } else {
        return ((static_cast<DefSFBitset>
            (z_ull - 1) & SFBitsetAux3D::k0SFBitsetTakeZRef_
            .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_
                .at(SFBitsetAux3D::kRefCurrent_))
            | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_
                .at(SFBitsetAux3D::kRefCurrent_)));
    }
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive z direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighboring node.
*/
inline DefSFBitset SFBitsetAux3D::FindZPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset sfbitest_tmp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((sfbitest_tmp & SFBitsetAux3D::k0SFBitsetTakeZRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code coincide with code at how many lower refinement levels.
* @param[in]  bitset_in   morton code of the current node.
* @return  number of coincident lower refinement levels.
*/
inline DefInt SFBitsetAux3D::SFBitsetCoincideLevel(const DefSFBitset& bitset_in) const {
    return static_cast<DefInt>(std::countr_zero(bitset_in.to_ullong()))/3;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_GRID_MORTON_AUX_HPP_
