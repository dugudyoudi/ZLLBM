//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.hpp
* @author Zhengliang Liu
* @date  2022-8-16
* @brief incline functions to manipulate morton code
*/

#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_MORTON_AUX_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_MORTON_AUX_H_
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to convert morton code to n level lower
* @param[in]  n_level  number of levels need to be shrinked.
* @param[in]  morton_in   input morton code.
* @return   morton code at lower levels.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToNLowerLevel (
    const DefSizet n_level,
    const DefSFBitset& morton_in)  const {
    return morton_in >> (2 * n_level);
}
/**
* @brief function to convert morton code to n level lower
* @param[in]  n_level  number of levels need to be streched.
* @param[in]  morton_in   input morton code.
* @return   morton code has been shrinked.
*/
inline DefSFBitset SFBitsetAux2D::SFBitsetToNHigherLevel(
    const DefSizet n_level,
    const DefSFBitset& morton_in) const {
    return morton_in << (2 * n_level);
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the negative x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux2D::FindXNeg(
    const DefSFBitset& bitset_in) const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefCurrent_)).to_ullong() - 1);
    return ((temp & SFBitsetAux2D::k0SFBitsetTakeXRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the positive x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux2D::FindXPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefOthers_)).to_ullong() + 1);
    return ((temp & SFBitsetAux2D::k0SFBitsetTakeXRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the negative y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux2D::FindYNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in & SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefCurrent_)).to_ullong() - 1);
    return ((temp & SFBitsetAux2D::k0SFBitsetTakeYRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
/**
* @brief function to find morton code (2D) of
*        the neighbour in the positive y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux2D::FindYPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux2D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux2D::kRefOthers_)).to_ullong() + 1);
    return ((temp & SFBitsetAux2D::k0SFBitsetTakeYRef_.
        at(SFBitsetAux2D::kRefCurrent_))
        | (bitset_in & SFBitsetAux2D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux2D::kRefCurrent_)));
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to convert morton code to n level lower
* @param[in]  n_level  number of levels need to be shrinked.
* @param[in]  morton_in   input morton code.
* @return   morton code at lower levels.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToNLowerLevel(
    const DefSizet n_level,
    const DefSFBitset& morton_in) const {
    return morton_in >> (3 * n_level);
}
/**
* @brief function to convert morton code to n level lower
* @param[in]  n_level  number of levels need to be streched.
* @param[in]  morton_in   input morton code.
* @return   morton code has been shrinked.
*/
inline DefSFBitset SFBitsetAux3D::SFBitsetToNHigherLevel(
    const DefSizet n_level,
    const DefSFBitset& morton_in)  const {
    return morton_in << (3 * n_level);
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the negative x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindXNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_)).to_ullong() - 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeXRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive x direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindXPos(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeXRef_.
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
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindYNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_)).to_ullong() - 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeYRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive y direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindYPos(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeYRef_.
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
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindZNeg(
    const DefSFBitset& bitset_in)  const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in & SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefCurrent_)).to_ullong() - 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeZRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
/**
* @brief function to find morton code (3D) of
*        the neighbour in the positive z direction.
* @param[in]  bitset_in   morton code of the current node.
* @return   morton code of the neighbouring node.
*/
inline DefSFBitset SFBitsetAux3D::FindZPos(
    const DefSFBitset& bitset_in) const {
    DefSFBitset temp = static_cast<DefSFBitset>
        ((bitset_in | SFBitsetAux3D::k0SFBitsetTakeZRef_.
            at(SFBitsetAux3D::kRefOthers_)).to_ullong() + 1);
    return ((temp & SFBitsetAux3D::k0SFBitsetTakeZRef_.
        at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeXRef_.
            at(SFBitsetAux3D::kRefCurrent_))
        | (bitset_in & SFBitsetAux3D::k0SFBitsetTakeYRef_.
            at(SFBitsetAux3D::kRefCurrent_)));
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_MORTON_AUX_H_
