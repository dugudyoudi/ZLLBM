//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file morton_aux.cpp
* @author Zhengliang Liu
* @brief functions used to manipulate morton codes.
* @date  2022-9-16
* @note  functions from other managers will be called.
*/
#include<array>
#include<string>
#include "auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
// static members
std::array<DefSFBitset, 2> SFBitsetAuxInterface::k0SFBitsetTakeXRef_, SFBitsetAuxInterface::k0SFBitsetTakeYRef_;
#ifndef  DEBUG_DISABLE_3D_FUNCTION
std::array<DefSFBitset, 2> SFBitsetAuxInterface::k0SFBitsetTakeZRef_;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to set bounds of 2D coordinates in bitset format.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
* @return  morton code.
*/
void SFBitsetAux2D::SFBitsetSetMinAndMaxBounds(
    const std::array<DefAmrIndexLUint, 2>& indices_min,
    const std::array<DefAmrIndexLUint, 2>& indices_max) {
    std::array<DefAmrIndexLUint, 2> indices_temp = { indices_min[kXIndex], 0 };
    SFBitsetMin_.at(kXIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { indices_max.at(kXIndex), 0 };
    SFBitsetMax_.at(kXIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { 0, indices_min[kYIndex] };
    SFBitsetMin_.at(kYIndex) = SFBitsetEncoding(indices_temp);
    indices_temp = { 0, indices_max.at(kYIndex) };
    SFBitsetMax_.at(kYIndex) = SFBitsetEncoding(indices_temp);
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to set bounds of 3D coordinates in bitset format.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
*/
void SFBitsetAux3D::SFBitsetSetMinAndMaxBounds(
    const std::array<DefAmrIndexLUint, 3>& indices_min,
    const std::array<DefAmrIndexLUint, 3>& indices_max) {
    SFBitsetMin_.at(kXIndex) =
        SFBitsetEncoding({ indices_min[kXIndex], 0, 0 });
    SFBitsetMax_.at(kXIndex) =
        SFBitsetEncoding({ indices_max[kXIndex], 0, 0 });
    SFBitsetMin_.at(kYIndex) =
        SFBitsetEncoding({ 0, indices_min[kYIndex], 0 });
    SFBitsetMax_.at(kYIndex) =
        SFBitsetEncoding({ 0, indices_max[kYIndex], 0 });
    SFBitsetMin_.at(kZIndex) =
        SFBitsetEncoding({ 0, 0, indices_min[kZIndex] });
    SFBitsetMax_.at(kZIndex) =
        SFBitsetEncoding({ 0, 0, indices_max[kZIndex] });
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to Identify if a node is not on the given
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
    std::array<bool, 2>* const ptr_bool_not_at_boundary_pos) const {
    // use reference sfbitset to take digits for x negative direction
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
*              sfbitset of the given node and its 3 neighbors.
* @note vec_sfbitsets[0]:(0, 0); vec_sfbitsets[1]:(+x, 0);
*       vec_sfbitsets[2]:(0, +y); vec_sfbitsets[3]:(+x, +y);
*/
void SFBitsetAux2D::SFBitsetFindCellNeighbors(
    const DefSFBitset& sfbitset_in,
    std::array<DefSFBitset, 4>* const ptr_sfbitsets)  const {
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
*              sfbitsets at the higher refinement level in a cell at the lower level.
*/
void SFBitsetAux2D::SFBitsetHigherLevelInACell(
    const DefAmrIndexUint level_diff,
    const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) const {
    SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet num = TwoPowerN(level_diff - 1);
    DefSFBitset sfbitset_x, sfbitset_y = sfbitset_corner;
    DefSizet yindex;
    for (DefSizet iy = 0; iy < num; ++iy) {
        yindex = iy * num;
        sfbitset_x = sfbitset_y;
        for (DefSizet ix = 0; ix < num; ++ix) {
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
* @param[out]  ptr_bitset_neighbors
*              sfbitset of the given node and its 8 neighbors.
*/
void SFBitsetAux2D::SFBitsetFindAllNeighbors(
    const DefSFBitset& sfbitset_center,
    std::array<DefSFBitset, 9>* const  ptr_bitset_neighbors) const {
    bool bool_x_gt0 = (sfbitset_center & k0SFBitsetTakeXRef_[kRefCurrent_]) != 0;
    bool bool_y_gt0 = (sfbitset_center & k0SFBitsetTakeYRef_[kRefCurrent_]) != 0;
    DefSFBitset sfbitset_temp0 = sfbitset_center, sfbitset_temp1 = sfbitset_center;

    ptr_bitset_neighbors->at(kNodeIndexX0Y0_) = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_temp0 = bool_x_gt0 ? FindXNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexXnY0_) = sfbitset_temp0;
    // node at (-x, -y, 0)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXnYn_) = sfbitset_temp1;
    // node at (-x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXnYp_) = sfbitset_temp1;
    // node at (+x, 0, 0)
    sfbitset_temp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexXpY0_) = sfbitset_temp0;
    // node at (+x, -y, 0)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXpYn_) = sfbitset_temp1;
    // node at (+x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXpYp_) = sfbitset_temp1;
    // node at (0, -y, 0)
    sfbitset_temp0 = bool_y_gt0 ? FindYNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0Yn_) = sfbitset_temp0;
    // node at (0, +y, 0)
    sfbitset_temp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0Yp_) = sfbitset_temp0;
}
/**
* @brief   function to find all neighboring nodes bonded by lower and upper limits.
* @param[in]  sfbitset_center sfbitset at the center.
* @param[in]  domain_min_n_level   lower bound.
* @param[in]  domain_max_n_level   upper bound.
* @param[out]  ptr_bitset_neighbors   space fill code for neighbors of the give node.
*/
void SFBitsetAux2D::SFBitsetFindAllBondedNeighborsVir(
    const DefSFBitset& bitset_in,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_bitset_neighbors) const {
    DefSFBitset sfbitset_temp0 = bitset_in;
    ptr_bitset_neighbors->clear();
    bool bool_not_x_min = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_min_n_level.at(kXIndex),
        bool_not_x_max = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_max_n_level.at(kXIndex),
        bool_not_y_min = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_min_n_level.at(kYIndex),
        bool_not_y_max = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_max_n_level.at(kYIndex);

    if (bool_not_x_min) {
        sfbitset_temp0 = FindXNeg(bitset_in);
        // node at (-x, 0, 0)
        ptr_bitset_neighbors->push_back(sfbitset_temp0);
        if (bool_not_y_min) {
            // node at (-x, -y, 0)
            ptr_bitset_neighbors->push_back(FindYNeg(sfbitset_temp0));
        }
        if (bool_not_y_max) {
            // node at (-x, +y, 0)
            ptr_bitset_neighbors->push_back(FindYPos(sfbitset_temp0));
        }
    }
    if (bool_not_x_max) {
        sfbitset_temp0 = FindXPos(bitset_in);
        // node at (+x, 0, 0)
        ptr_bitset_neighbors->push_back(sfbitset_temp0);
        if (bool_not_y_min) {
            // node at (+x, -y, 0)
            ptr_bitset_neighbors->push_back(FindYNeg(sfbitset_temp0));
        }
        if (bool_not_y_max) {
            // node at (+x, +y, 0)
            ptr_bitset_neighbors->push_back(FindYPos(sfbitset_temp0));
        }
    }
    if (bool_not_y_min) {
        // node at (0, -y, 0)
        ptr_bitset_neighbors->push_back(FindYNeg(bitset_in));
    }
    if (bool_not_y_max) {
        // node at (0, +y, 0)
        ptr_bitset_neighbors->push_back(FindYPos(bitset_in));
    }
}
/**
 * @brief function to search for the ghost layers near a given node based on min and max space fill codes.
 * @param[in] sfbitset_in space fill code of the given node
 * @param[in] region_length half length of the given region
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[out] ptr_sfbitset_nodes pointer to nodes in the given region.
 */
// example for a length 2 region, o is the input, * is negative x, x is positive x
// * * *  x  x
// * * *  x  x
// * * o  x  x
// * * *  x  x
// * * *  x  x
void SFBitsetAux2D::FindNodesInReginOfGivenLength(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length, const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level, std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    ptr_sfbitset_nodes->clear();
    // negative y direction
    for (DefAmrIndexLUint iy = 0; iy <= region_length; ++iy) {
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
         != domain_min_m1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                } else {
                    break;
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
        } else {
            break;
        }
        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
        sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
         != domain_max_p1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexLUint ix = 0; ix <= region_length; ++ix) {
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                } else {
                    break;
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
        } else {
            break;
        }
    }
}
/**
 * @brief function to find nodes in the region of a given length.
 * @param sfbitset_in the input node.
 * @param region_length the length of the region.
 * @param periodic_min  booleans indicating if the region is periodic in the minimum direction.
 * @param periodic_max booleans indicating if the region is periodic in the maximum direction.
 * @param domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid, otherwise may be undefined.
 */
// example for a length 2 region, o is the input, * is negative x, x is positive x
// * *  x  x
// * *  x  x
// * o  x  x
// * *  x  x
DefAmrIndexLUint SFBitsetAux2D::FindNodesInPeriodicReginOfGivenLength(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefAmrIndexLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, 0);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrIndexLUint vec_index_x, vec_index_y, index_min = region_length;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy - 1) * total_length + region_length;
        for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
            vec_index_x = vec_index_y - ix - 1;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > DefAmrIndexUint(ix + 1)) {
                        index_min = ix + 1;
                    }
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy - 1) * total_length + region_length;
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            } else {
                index_min = 0;
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                vec_index_x = vec_index_y + ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexUint(ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
            }
        }
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_min_n_level.at(kYIndex)) {
            if (periodic_min.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            } else {
                if (index_min > DefAmrIndexLUint(iy + 1)) {
                    index_min = iy + 1;
                }
                break;
            }
        }
        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    bool_not_y_max = true;
    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kYIndex)) {
        if (periodic_max.at(kYIndex)) {
            sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                |domain_min_n_level.at(kYIndex);
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
        } else {
            index_min = 0;
            bool_not_y_max = false;
        }
    }
    if (bool_not_y_max) {
        for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy) * total_length + region_length;
            for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > DefAmrIndexLUint(ix + 1)) {
                        index_min = ix + 1;
                    }
                    break;
                }
            }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy) * total_length + region_length;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    if (index_min > DefAmrIndexUint(iy + 1)) {
                        index_min = iy + 1;
                    }
                    break;
                }
            }
        }
    }
    return index_min;
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param sfbitset_in the input node.
 * @param region_length the length of the region.
 * @param periodic_min  booleans indicating if the region is periodic in the minimum direction.
 * @param periodic_max booleans indicating if the region is periodic in the maximum direction.
 * @param domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid, otherwise may be undefined.
 */
// example for a length 2 region, o is the input, * is nodes around it
// * * * * *
// * * * * *
// * * o * *
// * * * * *
// * * * * *
DefAmrIndexLUint SFBitsetAux2D::FindNodesInPeriodicReginNearby(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefAmrIndexLUint total_length = 2 * region_length + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, 0);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrIndexLUint vec_index_x, vec_index_y, index_min = region_length;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefAmrIndexLUint iy = 0; iy <= region_length; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy) * total_length + region_length;
        for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
            vec_index_x = vec_index_y - ix;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > DefAmrIndexUint(ix)) {
                        index_min = ix;
                    }
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy) * total_length + region_length;
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            } else {
                index_min = 0;
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                vec_index_x = vec_index_y + ix + 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexUint(ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
            }
        }
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_min_n_level.at(kYIndex)) {
            if (periodic_min.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            } else {
                if (index_min > DefAmrIndexLUint(iy)) {
                    index_min = iy;
                }
                break;
            }
        }
        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    bool_not_y_max = true;
    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kYIndex)) {
        if (periodic_max.at(kYIndex)) {
            sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                |domain_min_n_level.at(kYIndex);
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
        } else {
            index_min = 0;
            bool_not_y_max = false;
        }
    }
    if (bool_not_y_max) {
        for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy + 1) * total_length + region_length;
            for (DefAmrIndexLUint ix = 0; ix <= region_length; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > DefAmrIndexLUint(ix)) {
                        index_min = ix;
                    }
                    break;
                }
            }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy + 1) * total_length + region_length;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix + 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    if (index_min > DefAmrIndexUint(iy + 1)) {
                        index_min = iy + 1;
                    }
                    break;
                }
            }
        }
    }
    return index_min;
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 2
 */
void SFBitsetAux2D::GetMinM1AtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const {
    if (ptr_min_m1_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 2 in MpiManager::GetMinM1AtGivenLevel in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_min[kXIndex], 0}));
    ptr_min_m1_bitsets->at(kXIndex) = FindXNeg(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_min[kYIndex]}));
    ptr_min_m1_bitsets->at(kYIndex) = FindYNeg(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of maximum indices plus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[out] ptr_max_p1_bitsets a pointer to maximum indices plus 1
 * @throws ErrorType if the size of max_p1_bitsets is not 2
 */
void SFBitsetAux2D::GetMaxP1AtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const {
    if (ptr_max_p1_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 2 in MpiManager::GetMaxP1AtGivenLevel in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0}));
    ptr_max_p1_bitsets->at(kXIndex) = FindXPos(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex]}));
    ptr_max_p1_bitsets->at(kYIndex) = FindYPos(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of minimum indices at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_bitsets a pointer to minimum indices
 */
void SFBitsetAux2D::GetMinAtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_bitsets) const {
    if (ptr_min_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 2 in MpiManager::GetMinAtGivenLevel in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    ptr_min_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_min[kXIndex], 0}));
    ptr_min_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_min[kYIndex]}));
}
/**
 * @brief function to calculate spacing fill code of maximum indices at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[out] ptr_max_bitsets a pointer to maximum indices
 */
void SFBitsetAux2D::GetMaxAtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_bitsets) const {
    if (ptr_max_bitsets->size() != 2) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 2 in MpiManager::GetMaAtGivenLevel in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    ptr_max_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0}));
    ptr_max_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex]}));
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to Identify if a node is not on the given
*        cube boundary (3D).
* @param[in]  sfbitset_in  bitset of node need to be checked.
* @param[in]  sfbitset_min  bitset corresponding to the minimum coordinate
*                               in each direction.
* @param[in]  sfbitset_max  bitset corresponding to the maximum coordinate
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
    std::array<bool, 3>* const ptr_bool_not_at_boundary_pos) const {
    // use reference sfbitset to take digits for x negative direction
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
    const DefAmrIndexUint level_diff,
    const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) const {
    SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet num = TwoPowerN(level_diff - 1);
    DefSFBitset sfbitset_x, sfbitset_y, sfbitset_z = sfbitset_corner;
    DefSizet yindex, zindex;
    for (DefSizet iz = 0; iz < num; ++iz) {
        zindex = iz * num * num;
        for (DefSizet iy = 0; iy < num; ++iy) {
            yindex = zindex + iy * num;
            sfbitset_x = sfbitset_y;
            for (DefSizet ix = 0; ix < num; ++ix) {
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
*              sfbitset of the given node and its 7 neighbors.
* @note vec_sfbitsets[0]:(0, 0, 0); vec_sfbitsets[1]:(+x, 0, 0);
*       vec_sfbitsets[2]:(0, +y, 0); vec_sfbitsets[3]:(+x, +y, 0);
*       vec_sfbitsets[4]:(0, 0, +z);  vec_sfbitsets[5]:(+x, 0, +z);
*       vec_sfbitsets[6]:(0, +y, +z); vec_sfbitsets[7]: (+x, +y, +z).
*/
void SFBitsetAux3D::SFBitsetFindCellNeighbors(
    const DefSFBitset& sfbitset_in,
    std::array<DefSFBitset, 8>* const ptr_sfbitsets)  const {
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
* @brief   function to find all neighoring sfbitset (3D).
* @param[in]  sfbitset_center  sfbitset at the center.
* @param[out]  ptr_bitset_neighbors
*              sfbitset of the given node and its 26 neighbors.
*/
void SFBitsetAux3D::SFBitsetFindAllNeighbors(
    const DefSFBitset& sfbitset_center,
    std::array<DefSFBitset, 27>* const  ptr_bitset_neighbors) const {
    bool bool_x_gt0 = (sfbitset_center & k0SFBitsetTakeXRef_[kRefCurrent_]) != 0;
    bool bool_y_gt0 = (sfbitset_center & k0SFBitsetTakeYRef_[kRefCurrent_]) != 0;
    bool bool_z_gt0 = (sfbitset_center & k0SFBitsetTakeZRef_[kRefCurrent_]) != 0;
    DefSFBitset sfbitset_temp0 = sfbitset_center, sfbitset_temp1 = sfbitset_center,
     sfbitset_temp2 = sfbitset_center;

    ptr_bitset_neighbors->at(kNodeIndexX0Y0Z0_) = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_temp0 = bool_x_gt0 ? FindXNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Z0_) = sfbitset_temp0;
    // node at (-x, -y, 0)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZ0_) = sfbitset_temp1;
    // node at (-x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZ0_) = sfbitset_temp1;
    // node at (+x, 0, 0)
    sfbitset_temp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Z0_) = sfbitset_temp0;
    // node at (+x, -y, 0)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZ0_) = sfbitset_temp1;
    // node at (+x, +y, 0)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZ0_) = sfbitset_temp1;
    // node at (0, -y, 0)
    sfbitset_temp0 = bool_y_gt0 ? FindYNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZ0_) = sfbitset_temp0;
    // node at (0, +y, 0)
    sfbitset_temp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZ0_) = sfbitset_temp0;

    // node at (0, 0, -z)
    sfbitset_temp0 = bool_z_gt0 ? FindZNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0Y0Zn_) = sfbitset_temp0;
    // node at (-x, 0, -z)
    sfbitset_temp1 = bool_x_gt0 ? FindXNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Zn_) = sfbitset_temp1;
    // node at (-x, -y, -z)
    sfbitset_temp2 = bool_y_gt0 ? FindYNeg(sfbitset_temp1) : sfbitset_temp1;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZn_) = sfbitset_temp2;
    // node at (-x, +y, -z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZn_) = sfbitset_temp2;
    // node at (+x, 0, -z)
    sfbitset_temp1 = FindXPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Zn_) = sfbitset_temp1;
    // node at (+x, -y, -z)
    sfbitset_temp2 = bool_y_gt0 ? FindYNeg(sfbitset_temp1) : sfbitset_temp1;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZn_) = sfbitset_temp2;
    // node at (+x, +y, -z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZn_) = sfbitset_temp2;
    // node at (0, -y, -z)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZn_) = sfbitset_temp1;
    //  node at (0, +y, -z)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZn_) = sfbitset_temp1;
    // node at (0, 0, +z)
    sfbitset_temp0 = FindZPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0Y0Zp_) = sfbitset_temp0;
    // node at (-x, 0, +z)
    sfbitset_temp1 = bool_x_gt0 ? FindXNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Zp_) = sfbitset_temp1;
    // node at (-x, -y, +z)
    sfbitset_temp2 = bool_y_gt0 ? FindYNeg(sfbitset_temp1) : sfbitset_temp1;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZp_) = sfbitset_temp2;
    // node at (-x, +y, +z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZp_) = sfbitset_temp2;
    // node at (+x, 0, +z)
    sfbitset_temp1 = FindXPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Zp_) = sfbitset_temp1;
    // node at (+x, -y, +z)
    sfbitset_temp2 = bool_y_gt0 ? FindYNeg(sfbitset_temp1) : sfbitset_temp1;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZp_) = sfbitset_temp2;
    // node at (+x, +y, +z)
    sfbitset_temp2 = FindYPos(sfbitset_temp1);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZp_) = sfbitset_temp2;
    // node at (0, -y, +z)
    sfbitset_temp1 = bool_y_gt0 ? FindYNeg(sfbitset_temp0) : sfbitset_temp0;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZp_) = sfbitset_temp1;
    // node at (0, +y, +z)
    sfbitset_temp1 = FindYPos(sfbitset_temp0);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZp_) = sfbitset_temp1;
}
/**
* @brief   function to find all neighboring nodes bonded by lower and upper limits.
* @param[in]  sfbitset_center sfbitset at the center.
* @param[in]  domain_min_n_level   lower bound.
* @param[in]  domain_max_n_level   upper bound.
* @param[out]  ptr_bitset_neighbors   space fill code for neighbors of the give node.
*/
void SFBitsetAux3D::SFBitsetFindAllBondedNeighborsVir(
    const DefSFBitset& bitset_in,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_bitset_neighbors) const {
    DefSFBitset sfbitset_temp0 = bitset_in, sfbitset_temp1;
    ptr_bitset_neighbors->clear();
    bool bool_not_x_min = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_min_n_level.at(kXIndex),
        bool_not_x_max = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_max_n_level.at(kXIndex),
        bool_not_y_min = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_min_n_level.at(kYIndex),
        bool_not_y_max = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_max_n_level.at(kYIndex),
        bool_not_z_min = (bitset_in&k0SFBitsetTakeZRef_.at(kRefCurrent_)) != domain_min_n_level.at(kZIndex),
        bool_not_z_max = (bitset_in&k0SFBitsetTakeZRef_.at(kRefCurrent_)) != domain_max_n_level.at(kZIndex);

    if (bool_not_x_min) {
        sfbitset_temp0 = FindXNeg(bitset_in);
        // node at (-x, 0, 0)
        ptr_bitset_neighbors->push_back(sfbitset_temp0);
        if (bool_not_y_min) {
            // node at (-x, -y, 0)
            sfbitset_temp1 = FindYNeg(sfbitset_temp0);
            ptr_bitset_neighbors->push_back(sfbitset_temp1);
            if (bool_not_z_min) {
                // node at (-x, -y, -z)
                ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
            }
            if (bool_not_z_max) {
                // node at (-x, -y, +z)
                ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
            }
        }
        if (bool_not_y_max) {
            // node at (-x, +y, 0)
            sfbitset_temp1 = FindYPos(sfbitset_temp0);
            ptr_bitset_neighbors->push_back(sfbitset_temp1);
            if (bool_not_z_min) {
                // node at (-x, +y, -z)
                ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
            }
            if (bool_not_z_max) {
                // node at (-x, +y, +z)
                ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
            }
        }
        if (bool_not_z_min) {
                // node at (-x, 0, -z)
            ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp0));
        }
        if (bool_not_z_max) {
                // node at (-x, 0, +z)
            ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp0));
        }
    }
    if (bool_not_x_max) {
        sfbitset_temp0 = FindXPos(bitset_in);
        // node at (+x, 0, 0)
        ptr_bitset_neighbors->push_back(sfbitset_temp0);
        if (bool_not_y_min) {
            // node at (+x, -y, 0)
            sfbitset_temp1 = FindYNeg(sfbitset_temp0);
            ptr_bitset_neighbors->push_back(sfbitset_temp1);
            if (bool_not_z_min) {
                // node at (+x, -y, -z)
                ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
            }
            if (bool_not_z_max) {
                // node at (+-x, -y, +z)
                ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
            }
        }
        if (bool_not_y_max) {
            // node at (+x, +y, 0)
            sfbitset_temp1 = FindYPos(sfbitset_temp0);
            ptr_bitset_neighbors->push_back(sfbitset_temp1);
            if (bool_not_z_min) {
                // node at (+x, +y, -z)
                ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
            }
            if (bool_not_z_max) {
                // node at (+x, +y, +z)
                ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
            }
        }
        if (bool_not_z_min) {
                // node at (+x, 0, -z)
            ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp0));
        }
        if (bool_not_z_max) {
                // node at (+x, 0, +z)
            ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp0));
        }
    }
    if (bool_not_y_min) {
        // node at (0, -y, 0)
        sfbitset_temp1 = FindYNeg(bitset_in);
        ptr_bitset_neighbors->push_back(sfbitset_temp1);
        if (bool_not_z_min) {
            // node at (0, -y, -z)
            ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
        }
        if (bool_not_z_max) {
            // node at (0, -y, +z)
            ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
        }
    }
    if (bool_not_y_max) {
        // node at (0, +y, 0)
        sfbitset_temp1 = FindYPos(bitset_in);
        ptr_bitset_neighbors->push_back(sfbitset_temp1);
        if (bool_not_z_min) {
            // node at (0, +y, -z)
            ptr_bitset_neighbors->push_back(FindZNeg(sfbitset_temp1));
        }
        if (bool_not_z_max) {
            // node at (0, +y, +z)
            ptr_bitset_neighbors->push_back(FindZPos(sfbitset_temp1));
        }
    }
    if (bool_not_z_min) {
        // node at (0, 0, -z)
        ptr_bitset_neighbors->push_back(FindZNeg(bitset_in));
    }
    if (bool_not_z_max) {
        // node at (0, 0, +z)
        ptr_bitset_neighbors->push_back(FindZPos(bitset_in));
    }
}
/**
 * @brief function to search for the ghost layers near a given node based on min and max space fill codes.
 * @param[in] sfbitset_in space fill code of the given node
 * @param[in] region_length half length of the given region
 * @param[in] domain_min_m1_n_level minimum indicies of current refinement level minus 1.
 * @param[in] domain_max_p1_n_level maximum indicies of current refinement level plus 1.
 * @param[out] ptr_sfbitset_nodes pointer to nodes in the given region.
 */
void SFBitsetAux3D::FindNodesInReginOfGivenLength(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length, const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level, std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefSFBitset sfbitset_tmp_y, sfbitset_tmp_x, sfbitset_tmp_z = sfbitset_in;
    ptr_sfbitset_nodes->clear();
    // negative z direction
    for (DefAmrIndexUint iz = 0; iz <= region_length; ++iz) {
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
         != domain_min_m1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefAmrIndexUint iy = 0; iy <= region_length; ++iy) {
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefAmrIndexUint iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }
        } else {
            break;
        }
        sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    for (DefAmrIndexUint iz = 0; iz < region_length; ++iz) {
        sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
         != domain_max_p1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefAmrIndexUint iy = 0; iy <= region_length; ++iy) {
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            for (DefAmrIndexUint iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->push_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }
        } else {
            break;
        }
    }
}
/**
 * @brief function to search for the ghost layers near a given node based on min and max space fill codes.
 * @param[in] sfbitset_in space fill code of the given node.
 * @param[in] region_length half length of the given region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level minimum indicies of current refinement level.
 * @param[in] domain_max_n_level maximum indicies of current refinement level.
 * @param[out] ptr_sfbitset_nodes pointer to nodes in the given region.
 */
DefAmrIndexLUint SFBitsetAux3D::FindNodesInPeriodicReginOfGivenLength(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefAmrIndexLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, 0);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrIndexLUint vec_index_x, vec_index_y, vec_index_z, index_min = region_length;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefAmrIndexLUint iz = 0; iz < region_length; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length - iz - 1) * total_length  + region_length;
        for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
            for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexUint(ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kYIndex)) {
                if (periodic_min.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                } else {
                    if (index_min > DefAmrIndexLUint(iy + 1)) {
                        index_min = iy + 1;
                    }
                    break;
                }
            }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
        }
        // positive y direction
        sfbitset_tmp_y = sfbitset_tmp_z;
        bool_not_y_max = true;
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kYIndex)) {
            if (periodic_max.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy) * total_length + region_length;
                for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexLUint(ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy) * total_length + region_length;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                            } else {
                                if (index_min > DefAmrIndexUint(ix + 1)) {
                                    index_min = ix + 1;
                                }
                                break;
                            }
                        }
                    }
                }
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kYIndex)) {
                    if (periodic_max.at(kYIndex)) {
                        sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                    } else {
                        if (index_min > DefAmrIndexUint(iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                }
            }
        }
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kZIndex)) {
            if (periodic_min.at(kZIndex)) {
                sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kZIndex);
                sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            } else {
                if (index_min > DefAmrIndexLUint(iz + 1)) {
                    index_min = iz + 1;
                }
                break;
            }
        }
        sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    bool_not_z_max = true;
    if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kZIndex)) {
        if (periodic_max.at(kZIndex)) {
            sfbitset_tmp_z = (sfbitset_tmp_x&k0SFBitsetTakeZRef_.at(kRefOthers_))
                |domain_min_n_level.at(kZIndex);
            sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
        } else {
            index_min = 0;
            bool_not_z_max = false;
        }
    }
    if (bool_not_z_max) {
        for (DefAmrIndexLUint iz = 0; iz < region_length; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length + iz) * total_length + region_length;
            for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
                for (DefAmrIndexUint ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                            } else {
                                if (index_min > DefAmrIndexUint(ix + 1)) {
                                    index_min = ix + 1;
                                }
                                break;
                            }
                        }
                    }
                }
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kYIndex)) {
                    if (periodic_min.at(kYIndex)) {
                        sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    } else {
                        if (index_min > DefAmrIndexLUint(iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            bool_not_y_max = true;
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy) * total_length + region_length;
                    for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                        vec_index_x = vec_index_y - ix - 1;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexLUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy) * total_length + region_length;
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            index_min = 0;
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                            vec_index_x = vec_index_y + ix;
                            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                                == domain_max_n_level.at(kXIndex)) {
                                if (periodic_max.at(kXIndex)) {
                                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                                } else {
                                    if (index_min > DefAmrIndexUint(ix + 1)) {
                                        index_min = ix + 1;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kYIndex)) {
                        if (periodic_max.at(kYIndex)) {
                            sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                        } else {
                            if (index_min > DefAmrIndexUint(iy + 1)) {
                                index_min = iy + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kZIndex)) {
                if (periodic_max.at(kZIndex)) {
                    sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kZIndex);
                    sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
                } else {
                    if (index_min > DefAmrIndexUint(iz + 1)) {
                        index_min = iz + 1;
                    }
                    break;
                }
            }
        }
    }
    return index_min;
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param sfbitset_in the input node.
 * @param region_length the length of the region.
 * @param periodic_min  booleans indicating if the region is periodic in the minimum direction.
 * @param periodic_max booleans indicating if the region is periodic in the maximum direction.
 * @param domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid, otherwise may be undefined.
 */
DefAmrIndexLUint SFBitsetAux3D::FindNodesInPeriodicReginNearby(const DefSFBitset& sfbitset_in,
    const DefAmrIndexLUint region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
    DefAmrIndexLUint total_length = 2 * region_length + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, 0);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrIndexLUint vec_index_x, vec_index_y, vec_index_z, index_min = region_length;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefAmrIndexLUint iz = 0; iz <= region_length; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length - iz) * total_length  + region_length;
        for (DefAmrIndexLUint iy = 0; iy <= region_length; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length + region_length;
            for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexUint(ix)) {
                            index_min = ix;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length + region_length;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix + 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kYIndex)) {
                if (periodic_min.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                } else {
                    if (index_min > DefAmrIndexLUint(iy)) {
                        index_min = iy;
                    }
                    break;
                }
            }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
        }
        // positive y direction
        sfbitset_tmp_y = sfbitset_tmp_z;
        bool_not_y_max = true;
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kYIndex)) {
            if (periodic_max.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                for (DefAmrIndexLUint ix = 0; ix <= region_length; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefAmrIndexLUint(ix)) {
                            index_min = ix;
                        }
                        break;
                    }
                }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix + 1;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                            } else {
                                if (index_min > DefAmrIndexUint(ix + 1)) {
                                    index_min = ix + 1;
                                }
                                break;
                            }
                        }
                    }
                }
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kYIndex)) {
                    if (periodic_max.at(kYIndex)) {
                        sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                    } else {
                        if (index_min > DefAmrIndexUint(iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                }
            }
        }
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kZIndex)) {
            if (periodic_min.at(kZIndex)) {
                sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kZIndex);
                sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            } else {
                if (index_min > DefAmrIndexLUint(iz)) {
                    index_min = iz;
                }
                break;
            }
        }
        sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    bool_not_z_max = true;
    if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kZIndex)) {
        if (periodic_max.at(kZIndex)) {
            sfbitset_tmp_z = (sfbitset_tmp_x&k0SFBitsetTakeZRef_.at(kRefOthers_))
                |domain_min_n_level.at(kZIndex);
            sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
        } else {
            index_min = 0;
            bool_not_z_max = false;
        }
    }
    if (bool_not_z_max) {
        for (DefAmrIndexLUint iz = 0; iz < region_length; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length + iz + 1) * total_length + region_length;
            for (DefAmrIndexLUint iy = 0; iy <= region_length; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length + region_length;
                for (DefAmrIndexUint ix = 0; ix <= region_length; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexUint(ix)) {
                                index_min = ix;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length + region_length;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix + 1;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                            } else {
                                if (index_min > DefAmrIndexUint(ix + 1)) {
                                    index_min = ix + 1;
                                }
                                break;
                            }
                        }
                    }
                }
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kYIndex)) {
                    if (periodic_min.at(kYIndex)) {
                        sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    } else {
                        if (index_min > DefAmrIndexLUint(iy)) {
                            index_min = iy;
                        }
                        break;
                    }
                }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            bool_not_y_max = true;
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_x&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefAmrIndexLUint iy = 0; iy < region_length; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                    for (DefAmrIndexLUint ix = 0; ix <= region_length; ++ix) {
                        vec_index_x = vec_index_y - ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > DefAmrIndexLUint(ix)) {
                                index_min = ix;
                            }
                            break;
                        }
                    }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            index_min = 0;
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefAmrIndexLUint ix = 0; ix < region_length; ++ix) {
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                            vec_index_x = vec_index_y + ix + 1;
                            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                                == domain_max_n_level.at(kXIndex)) {
                                if (periodic_max.at(kXIndex)) {
                                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                                } else {
                                    if (index_min > DefAmrIndexUint(ix + 1)) {
                                        index_min = ix + 1;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kYIndex)) {
                        if (periodic_max.at(kYIndex)) {
                            sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                        } else {
                            if (index_min > DefAmrIndexUint(iy + 1)) {
                                index_min = iy + 1;
                            }
                            break;
                        }
                    }
                }
            }
            if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kZIndex)) {
                if (periodic_max.at(kZIndex)) {
                    sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kZIndex);
                    sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
                } else {
                    if (index_min > DefAmrIndexUint(iz + 1)) {
                        index_min = iz + 1;
                    }
                    break;
                }
            }
        }
    }
    return index_min;
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 3
 */
void SFBitsetAux3D::GetMinM1AtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const {
    if (ptr_min_m1_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 3 in MpiManager::GetMinM1AtGivenLevel3D in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_min[kXIndex], 0, 0}));
    ptr_min_m1_bitsets->at(kXIndex) = FindXNeg(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0,  indices_min[kYIndex], 0}));
    ptr_min_m1_bitsets->at(kYIndex) = FindYNeg(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, 0, indices_min[kZIndex]}));
    ptr_min_m1_bitsets->at(kZIndex) = FindZNeg(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of maximum indices plus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[out] ptr_max_p1_bitsets a pointer to maximum indices plus 1
 * @throws ErrorType if the size of max_p1_bitsets is not 3
 */
void SFBitsetAux3D::GetMaxP1AtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const {
    if (ptr_max_p1_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 3 in MpiManager::GetMaxP1AtGivenLevel3D in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFBitset bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0, 0}));
    ptr_max_p1_bitsets->at(kXIndex) = FindXPos(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex], 0}));
    ptr_max_p1_bitsets->at(kYIndex) = FindYPos(bitset_tmp);
    bitset_tmp = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, 0, indices_max[kZIndex]}));
    ptr_max_p1_bitsets->at(kZIndex) = FindZPos(bitset_tmp);
}
/**
 * @brief function to calculate spacing fill code of minimum indices at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_bitsets a pointer to minimum indices
 */
void SFBitsetAux3D::GetMinAtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_bitsets) const {
    if (ptr_min_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_min_m1_bitsets should be 3 in MpiManager::GetMinAtGivenLevel3D in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    ptr_min_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_min[kXIndex], 0, 0}));
    ptr_min_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0,  indices_min[kYIndex], 0}));
    ptr_min_bitsets->at(kZIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, 0, indices_min[kZIndex]}));
}
/**
 * @brief function to calculate spacing fill code of maximum indices at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[out] ptr_max_bitsets a pointer to maximum indices
 */
void SFBitsetAux3D::GetMaxAtGivenLevel(const DefAmrIndexUint i_level,
    std::vector<DefAmrIndexLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_bitsets) const {
    if (ptr_max_bitsets->size() != 3) {
        LogManager::LogError("size of ptr_max_p1_bitsets should be 3 in MpiManager::GetMaxAtGivenLevel3D in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    ptr_max_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0, 0}));
    ptr_max_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex], 0}));
    ptr_max_bitsets->at(kZIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, 0, indices_max[kZIndex]}));
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject