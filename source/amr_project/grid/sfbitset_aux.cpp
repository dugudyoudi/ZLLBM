//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.cpp
* @author Zhengliang Liu
* @brief functions used to manipulate morton codes.
* @date  2022-9-16
* @note  functions from other managers will be called.
*/
#include<array>
#include "auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
namespace grid {
// static members
#ifndef  DEBUG_DISABLE_2D_FUNCTION
std::array<DefSFBitset, 2> SFBitsetAuxInterface::k0SFBitsetTakeXRef_,
SFBitsetAuxInterface::k0SFBitsetTakeYRef_;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
std::array<DefSFBitset, 2> SFBitsetAuxInterface::k0SFBitsetTakeZRef_;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to set bounds of 2D coordinates in bitset format.
* @param[in]  max_level      maximum refinement level.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
* @return  morton code.
*/
void SFBitsetAux2D::SFBitsetMinAndMaxCoordinates(
    const DefSizet max_level,
    const std::array<DefLUint, 2>& indices_min,
    const std::array<DefLUint, 2>& indices_max) {
    std::array<DefLUint, 2> indices_temp = { indices_min[kXIndex], 0 };
    k0SFBitsetMin.at(kXIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { indices_max.at(kXIndex), 0 };
    k0SFBitsetMax.at(kXIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { 0, indices_min[kYIndex] };
    k0SFBitsetMin.at(kYIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { 0, indices_max.at(kYIndex) };
    k0SFBitsetMax.at(kYIndex) = SFBitsetEncoding(indices_temp);
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS  
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to set bounds of 3D coordinates in bitset format.
* @param[in]  max_level      maximum refinement level.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
* @return  morton code.
*/
void SFBitsetAux3D::SFBitsetMinAndMaxCoordinates(
    const DefSizet max_level,
    const std::array<DefLUint, 3>& indices_min,
    const std::array<DefLUint, 3>& indices_max) {
        k0SFBitsetMin.at(kXIndex) =
            SFBitsetEncoding({ indices_min[kXIndex], 0, 0 });
        k0SFBitsetMax.at(kXIndex) =
            SFBitsetEncoding({ indices_max[kXIndex], 0, 0 });
        k0SFBitsetMin.at(kYIndex) =
            SFBitsetEncoding({ 0, indices_min[kYIndex], 0 });
        k0SFBitsetMax.at(kYIndex) =
            SFBitsetEncoding({ 0, indices_max[kYIndex], 0 });
        k0SFBitsetMin.at(kZIndex) =
            SFBitsetEncoding({ 0, 0, indices_min[kZIndex] });
        k0SFBitsetMax.at(kZIndex) =
            SFBitsetEncoding({ 0, 0, indices_max[kZIndex] });
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to indentify if a node is not on the given
*        rectangle boundary (2D).
* @param[in]  sfbitset_in  bitset of node need to be checked.
* @param[in]  sfbitset_min  bitset corresponding to the minimum coordinate
*                               in each direction.
* @param[in]  sfbitset_max  bitset corresponding to the maximim coordinate
*                               in each direction. 
* @param[out]  ptr_bool_not_at_boundary_neg
*                identifier of node on the boundary in negative directions.
* @param[out]  ptr_bool_not_at_boundary_pos
*                identifier of node on the boundary in positive directions.
*/
void SFBitsetAux2D::SFBitsetNotOnDomainBoundary(
    const DefSFBitset& sfbitset_in,
    const std::array<DefSFBitset, 2>& sfbitset_min,
    const std::array<DefSFBitset, 2>& sfbitset_max,
    std::array<bool, 2>* const ptr_bool_not_at_boundary_neg,
    std::array<bool, 2>* const ptr_bool_not_at_boundary_pos) {
    // use refernece sfbitset to take digits for x negative direction
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_min[kXIndex]) {
        ptr_bool_not_at_boundary_neg->at(kXIndex) = false;
    } else {
        ptr_bool_not_at_boundary_neg->at(kXIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_max[kXIndex]) {
        ptr_bool_not_at_boundary_pos->at(kXIndex) = false;
    } else {
        ptr_bool_not_at_boundary_pos->at(kXIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_min[kYIndex]) {
        ptr_bool_not_at_boundary_neg->at(kYIndex) = false;
    } else {
        ptr_bool_not_at_boundary_neg->at(kYIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_max[kYIndex]) {
        ptr_bool_not_at_boundary_pos->at(kYIndex) = false;
    } else {
        ptr_bool_not_at_boundary_pos->at(kYIndex) = true;
    }
}
/**
* @brief   function to calculate sfbitset constructing a cell (2D).
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_sfbitsets
*              sfbitset of the given node and its 3 neighbours.
* @note vec_sfbitsets[0]:(0, 0); vec_sfbitsets[1]:(+x, 0);
*       vec_sfbitsets[2]:(0, +y); vec_sfbitsets[3]:(+x, +y);
*/
void SFBitsetAux2D::SFBitsetFindCellNeighbours(
    const DefSFBitset& sfbitset_in,
    std::array<DefSFBitset, 4>* const ptr_sfbitsets) {
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    // (+x, +y)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    // (0, +y)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
}
/**
* @brief   function to find all sfbitset at n higher level (2D).
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_vec_sfbitsets_higher_level
*              sfbitsets at the higher refinment level in a cell
*              at the lower level.
*/
void SFBitsetAux2D::SFBitsetHigherLevelInACell(
    const DefSizet level_diff,
    const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) {
    SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet num = TwoPowerN(level_diff - 1);
    DefSFBitset sfbitset_x, sfbitset_y = sfbitset_corner;
    DefSizet yindex;
    for (DefSizet iy = 0; iy < num; ++iy) {
        yindex = iy * num;
        sfbitset_x = sfbitset_y;
        for (DefUint ix = 0; ix < num; ++ix) {
            ptr_vec_sfbitsets_higher_level->at(yindex + ix)
                = sfbitset_x;
            sfbitset_x = FindXPos(sfbitset_x);
        }
        sfbitset_y = FindYPos(sfbitset_y);
    }
}
/**
* @brief   function to find all neighouring sfbitset (2D).
* @param[in]  sfbitset_center sfbitset at the center.
* @param[out]  ptr_bitset_neighbours
*              sfbitset of the given node and its 8 neighbours.
*/
void SFBitsetAux2D::SFBitsetFindAllNeighbours(
    const DefSFBitset& sfbitset_center,
    std::array<DefSFBitset, 9>* const  ptr_bitset_neighbours) {

    DefSFBitset sfbitset_temp0, sfbitset_temp1;

    ptr_bitset_neighbours->at(kNodeIndexX0Y0_)
        = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_temp0 = FindXNeg(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexXnY0_)
        = sfbitset_temp0;
    // node at (-x, -y, 0)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnYn_)
        = sfbitset_temp1;
    // node at (-x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnYp_)
        = sfbitset_temp1;
    // node at (+x, 0, 0)
    sfbitset_temp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexXpY0_)
        = sfbitset_temp0;
    // node at (+x, -y, 0)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpYn_)
        = sfbitset_temp1;
    // node at (+x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpYp_)
        = sfbitset_temp1;
    // node at (0, -y, 0)
    sfbitset_temp0 = FindYNeg(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0Yn_)
        = sfbitset_temp0;
    // node at (0, +y, 0)
    sfbitset_temp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0Yp_)
        = sfbitset_temp0;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to indentify if a node is not on the given
*        cube boundary (3D).
* @param[in]  sfbitset_in  bitset of node need to be checked.
* @param[in]  sfbitset_min  bitset corresponding to the minimum coordinate
*                               in each direction.
* @param[in]  sfbitset_max  bitset corresponding to the maximim coordinate
*                               in each direction.
* @param[out]  ptr_bool_not_at_boundary_neg
*                identifier of node on the boundary in negative directions.
* @param[out]  ptr_bool_not_at_boundary_pos
*                identifier of node on the boundary in positive directions.
*/
void SFBitsetAux3D::SFBitsetNotOnDomainBoundary(
    const DefSFBitset& sfbitset_in,
    const std::array<DefSFBitset, 3>& sfbitset_min,
    const std::array<DefSFBitset, 3>& sfbitset_max,
    std::array<bool, 3>* const ptr_bool_not_at_boundary_neg,
    std::array<bool, 3>* const ptr_bool_not_at_boundary_pos) {
    // use refernece sfbitset to take digits for x negative direction
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_min[kXIndex]) {
        ptr_bool_not_at_boundary_neg->at(kXIndex) = false;
    } else {
        ptr_bool_not_at_boundary_neg->at(kXIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_max[kXIndex]) {
        ptr_bool_not_at_boundary_pos->at(kXIndex) = false;
    } else {
        ptr_bool_not_at_boundary_pos->at(kXIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_min[kYIndex]) {
        ptr_bool_not_at_boundary_neg->at(kYIndex) = false;
    } else {
        ptr_bool_not_at_boundary_neg->at(kYIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_max[kYIndex]) {
        ptr_bool_not_at_boundary_pos->at(kYIndex) = false;
    } else {
        ptr_bool_not_at_boundary_pos->at(kYIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeZRef_[kRefCurrent_])
        == sfbitset_min[kZIndex]) {
        ptr_bool_not_at_boundary_neg->at(kZIndex) = false;
    } else {
        ptr_bool_not_at_boundary_neg->at(kZIndex) = true;
    }
    if ((sfbitset_in & k0SFBitsetTakeZRef_[kRefCurrent_])
        == sfbitset_max[kZIndex]) {
        ptr_bool_not_at_boundary_pos->at(kZIndex) = false;
    } else {
        ptr_bool_not_at_boundary_pos->at(kZIndex) = true;
    }
}
/**
* @brief   function to find all sfbitset at n higher level (3D).
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_vec_sfbitsets_higher_level
*              sfbitsets at the higher refinment level in a cell
*              at the lower level.
*/
void SFBitsetAux3D::SFBitsetHigherLevelInACell(
    const DefSizet level_diff,
    const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) {
    SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet num = TwoPowerN(level_diff - 1);
    DefSFBitset sfbitset_x, sfbitset_y, sfbitset_z = sfbitset_corner;
    DefSizet yindex, zindex;
    for (DefSizet iz = 0; iz < num; ++iz) {
        zindex = iz * num * num;
        for (DefSizet iy = 0; iy < num; ++iy) {
            yindex = zindex + iy * num;
            sfbitset_x = sfbitset_y;
            for (DefUint ix = 0; ix < num; ++ix) {
                ptr_vec_sfbitsets_higher_level->at(yindex + ix)
                    = sfbitset_x;
                sfbitset_x = FindXPos(sfbitset_x);
            }
            sfbitset_y = FindYPos(sfbitset_y);
        }
        sfbitset_z = FindZPos(sfbitset_z);
    }
}
/**
* @brief   function to calculate sfbitset constructing a cell (3D).
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_sfbitsets
*              sfbitset of the given node and its 7 neighbours.
* @note vec_sfbitsets[0]:(0, 0, 0); vec_sfbitsets[1]:(+x, 0, 0);
*       vec_sfbitsets[2]:(0, +y, 0); vec_sfbitsets[3]:(+x, +y, 0);
*       vec_sfbitsets[4]:(0, 0, +z);  vec_sfbitsets[5]:(+x, 0, +z);
*       vec_sfbitsets[6]:(0, +y, +z); vec_sfbitsets[7]: (+x, +y, +z).
*/
void SFBitsetAux3D::SFBitsetFindCellNeighbours(
    const DefSFBitset& sfbitset_in,
    std::array<DefSFBitset, 8>* const ptr_sfbitsets) {
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    // (+x, +y, 0)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    // (0, +y, 0)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
    // (0, 0, +z)
    ptr_sfbitsets->at(4) = FindZPos(sfbitset_in);
    // (+x, 0, +z)
    ptr_sfbitsets->at(5) = FindXPos(ptr_sfbitsets->at(4));
    // (+x, +y, +z)
    ptr_sfbitsets->at(7) = FindYPos(ptr_sfbitsets->at(5));
    // (0, +y, +z)
    ptr_sfbitsets->at(6) = FindYPos(ptr_sfbitsets->at(4));
}
/**
* @brief   function to find all neighouring sfbitset (3D).
* @param[in]  sfbitset_center  sfbitset at the center.
* @param[out]  ptr_bitset_neighbours
*              sfbitset of the given node and its 26 neighbours.
*/
void SFBitsetAux3D::SFBitsetFindAllNeighbours(
    const DefSFBitset& sfbitset_center,
    std::array<DefSFBitset,27>* const  ptr_bitset_neighbours) {

    DefSFBitset sfbitset_temp0, sfbitset_temp1, sfbitset_temp2;

    ptr_bitset_neighbours->at(kNodeIndexX0Y0Z0_)
        = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_temp0 = FindXNeg(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexXnY0Z0_)
        = sfbitset_temp0;
    // node at (-x, -y, 0)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnYnZ0_)
        = sfbitset_temp1;
    // node at (-x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnYpZ0_)
        = sfbitset_temp1;
    // node at (+x, 0, 0)
    sfbitset_temp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexXpY0Z0_)
        = sfbitset_temp0;
    // node at (+x, -y, 0)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpYnZ0_)
        = sfbitset_temp1;
    // node at (+x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpYpZ0_)
        = sfbitset_temp1;
    // node at (0, -y, 0)
    sfbitset_temp0 = FindYNeg(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0YnZ0_)
        = sfbitset_temp0;
    // node at (0, +y, 0)
    sfbitset_temp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0YpZ0_)
        = sfbitset_temp0;

    // node at (0, 0, -z)
    sfbitset_temp0 = FindZNeg(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0Y0Zn_)
        = sfbitset_temp0;
    // node at (-x, 0, -z)
    sfbitset_temp1 = FindXNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnY0Zn_)
        = sfbitset_temp1;
    // node at (-x, -y, -z)
    sfbitset_temp2 = FindYNeg(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXnYnZn_)
        = sfbitset_temp2;
    // node at (-x, +y, -z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXnYpZn_)
        = sfbitset_temp2;
    // node at (+x, 0, -z)
    sfbitset_temp1 = FindXPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpY0Zn_)
        = sfbitset_temp1;
    // node at (+x, -y, -z)
    sfbitset_temp2 = FindYNeg(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXpYnZn_)
        = sfbitset_temp2;
    // node at (+x, +y, -z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXpYpZn_)
        = sfbitset_temp2;
    // node at (0, -y, -z)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexX0YnZn_)
        = sfbitset_temp1;
    //  node at (0, +y, -z)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexX0YpZn_)
        = sfbitset_temp1;
    // node at (0, 0, +z)
    sfbitset_temp0 = FindZPos(sfbitset_center);
    ptr_bitset_neighbours->at(kNodeIndexX0Y0Zp_)
        = sfbitset_temp0;
    // node at (-x, 0, +z)
    sfbitset_temp1 = FindXNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXnY0Zp_)
        = sfbitset_temp1;
    // node at (-x, -y, +z)
    sfbitset_temp2 = FindYNeg(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXnYnZp_)
        = sfbitset_temp2;
    // node at (-x, +y, +z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXnYpZp_)
        = sfbitset_temp2;
    // node at (+x, 0, +z)
    sfbitset_temp1 = FindXPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexXpY0Zp_)
        = sfbitset_temp1;
    // node at (+x, -y, +z)
    sfbitset_temp2 = FindYNeg(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXpYnZp_)
        = sfbitset_temp2;
    // node at (+x, +y, +z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbours->at(kNodeIndexXpYpZp_)
        = sfbitset_temp2;
    // node at (0, -y, +z)
    sfbitset_temp1 = FindYNeg(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexX0YnZp_)
        = sfbitset_temp1;
    //  node at (0, +y, +z)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbours->at(kNodeIndexX0YpZp_)
        = sfbitset_temp1;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namsapce grid
}  // end namespace amrproject
}  // end namespace rootproject