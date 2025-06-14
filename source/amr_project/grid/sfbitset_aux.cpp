//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
#include <utility>
#include "./auxiliary_inline_func.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
void SFBitsetAux2D::SetSpaceBackground(const std::vector<DefReal>& space_background) {
    if (space_background.size() != 2) {
        LogManager::LogError("dimension of space background should be 2");
    }
    k0SpaceBackground_ = space_background;
}
/**
* @brief function to set bounds of 2D coordinates in bitset format.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
* @return  morton code.
*/
void SFBitsetAux2D::SFBitsetSetMinAndMaxBounds(
    const std::array<DefAmrLUint, 2>& indices_min,
    const std::array<DefAmrLUint, 2>& indices_max) {
    std::array<DefAmrLUint, 2> indices_tmp = { indices_min[kXIndex], 0 };
    SFBitsetMin_.at(kXIndex) = SFBitsetEncoding(indices_tmp);
    indices_tmp = { indices_max.at(kXIndex), 0 };
    SFBitsetMax_.at(kXIndex) = SFBitsetEncoding(indices_tmp);
    indices_tmp = { 0, indices_min[kYIndex] };
    SFBitsetMin_.at(kYIndex) = SFBitsetEncoding(indices_tmp);
    indices_tmp = { 0, indices_max.at(kYIndex) };
    SFBitsetMax_.at(kYIndex) = SFBitsetEncoding(indices_tmp);
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief function to set bounds of 3D coordinates in bitset format.
* @param[in]  indices_min      minimum indices.
* @param[in]  indices_max      maximum indices.
*/
void SFBitsetAux3D::SFBitsetSetMinAndMaxBounds(
    const std::array<DefAmrLUint, 3>& indices_min,
    const std::array<DefAmrLUint, 3>& indices_max) {
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
* @param[in]  level_diff  different between target and given levels of nodes.
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_vec_sfbitsets_higher_level
*              sfbitsets at the higher refinement level in a cell at the lower level.
*/
void SFBitsetAux2D::SFBitsetHigherLevelInACell(
    const DefInt level_diff, const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) const {
    DefSizet num = TwoPowerN(level_diff);
    ptr_vec_sfbitsets_higher_level->resize(num*num);
    DefSFBitset sfbitset_x, sfbitset_y = SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet yindex;
    for (DefSizet iy = 0; iy < num; ++iy) {
        yindex = iy * num;
        sfbitset_x = sfbitset_y;
        for (DefSizet ix = 0; ix < num; ++ix) {
            ptr_vec_sfbitsets_higher_level->at(yindex + ix) = sfbitset_x;
            sfbitset_x = FindXPos(sfbitset_x);
        }
        sfbitset_y = FindYPos(sfbitset_y);
    }
}
/**
* @brief   function to find all neighbouring sfbitset (2D).
* @param[in]  sfbitset_center sfbitset at the center.
* @param[out]  ptr_bitset_neighbors
*              sfbitset of the given node and its 8 neighbors.
*/
void SFBitsetAux2D::SFBitsetFindAllNeighbors(
    const DefSFBitset& sfbitset_center,
    std::array<DefSFBitset, 9>* const  ptr_bitset_neighbors) const {
    bool bool_x_gt0 = (sfbitset_center & k0SFBitsetTakeXRef_[kRefCurrent_]) != 0;
    bool bool_y_gt0 = (sfbitset_center & k0SFBitsetTakeYRef_[kRefCurrent_]) != 0;
    DefSFBitset sfbitset_tmp0 = sfbitset_center, sfbitset_tmp1 = sfbitset_center;

    ptr_bitset_neighbors->at(kNodeIndexX0Y0_) = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_tmp0 = bool_x_gt0 ? FindXNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexXnY0_) = sfbitset_tmp0;
    // node at (-x, -y, 0)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXnYn_) = sfbitset_tmp1;
    // node at (-x, +y, 0)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXnYp_) = sfbitset_tmp1;
    // node at (+x, 0, 0)
    sfbitset_tmp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexXpY0_) = sfbitset_tmp0;
    // node at (+x, -y, 0)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXpYn_) = sfbitset_tmp1;
    // node at (+x, +y, 0)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXpYp_) = sfbitset_tmp1;
    // node at (0, -y, 0)
    sfbitset_tmp0 = bool_y_gt0 ? FindYNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0Yn_) = sfbitset_tmp0;
    // node at (0, +y, 0)
    sfbitset_tmp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0Yp_) = sfbitset_tmp0;
}
/**
* @brief   function to find all neighboring nodes bonded by lower and upper limits.
* @param[in]  sfbitset_center sfbitset at the center.
* @param[in] bool_periodic_min booleans indicating if the boundary is periodic at maximum domain boundaries.
* @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
* @param[in]  domain_min_n_level   lower bound.
* @param[in]  domain_max_n_level   upper bound.
* @param[out]  ptr_bitset_neighbors   space fill code for neighbors of the give node.
*/
void SFBitsetAux2D::SFBitsetFindAllBondedNeighborsVir(const DefSFBitset& bitset_in,
    const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level, const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_bitset_neighbors) const {
    DefSFBitset sfbitset_tmp0 = bitset_in;
    ptr_bitset_neighbors->clear();
    bool bool_not_x_min = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_min_n_level.at(kXIndex),
        bool_not_x_max = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_max_n_level.at(kXIndex),
        bool_not_y_min = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_min_n_level.at(kYIndex),
        bool_not_y_max = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_max_n_level.at(kYIndex);

    if (bool_not_x_min) {
        sfbitset_tmp0 = FindXNeg(bitset_in);
        // node at (-x, 0, 0)
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    } else if (bool_periodic_min.at(kXIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))|domain_max_n_level.at(kXIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_y_min) {
        // node at (-x, -y, 0)
        ptr_bitset_neighbors->emplace_back(FindYNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex));
    }
    if (bool_not_y_max) {
        // node at (-x, +y, 0)
        ptr_bitset_neighbors->emplace_back(FindYPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex));
    }

    if (bool_not_x_max) {
        sfbitset_tmp0 = FindXPos(bitset_in);
        // node at (+x, 0, 0)
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    } else if (bool_periodic_max.at(kXIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))|domain_min_n_level.at(kXIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_y_min) {
        // node at (+x, -y, 0)
        ptr_bitset_neighbors->emplace_back(FindYNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex));
    }
    if (bool_not_y_max) {
        // node at (+x, +y, 0)
        ptr_bitset_neighbors->emplace_back(FindYPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex));
    }

    if (bool_not_y_min) {
        // node at (0, -y, 0)
        ptr_bitset_neighbors->emplace_back(FindYNeg(bitset_in));
    } else if (bool_periodic_min.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (bitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex));
    }
    if (bool_not_y_max) {
        // node at (0, +y, 0)
        ptr_bitset_neighbors->emplace_back(FindYPos(bitset_in));
    } else if (bool_periodic_max.at(kYIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (bitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex));
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
    const DefInt region_length, const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level, std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_m1_n_level.size() != 2 || domain_max_p1_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_m1_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_p1_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    ptr_sfbitset_nodes->clear();
    // negative y direction
    for (DefInt iy = 0; iy <= region_length; ++iy) {
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
         != domain_min_m1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefInt ix = 0; ix <= region_length; ++ix) {
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                } else {
                    break;
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
    for (DefInt iy = 0; iy < region_length; ++iy) {
        sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
         != domain_max_p1_n_level.at(kYIndex)) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefInt ix = 0; ix <= region_length; ++ix) {
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                } else {
                    break;
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kXIndex)) {
                    ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length the length of the region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid.
 */
// example for a length 2 region, o is the input, * is negative x, x is positive x
// * *  x  x
// * *  x  x
// * o  x  x
// * *  x  x
DefInt SFBitsetAux2D::FindNodesInPeriodicRegionCorner(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy < region_length; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy - 1) * total_length + region_length;
        for (DefInt ix = 0; ix < region_length; ++ix) {
            vec_index_x = vec_index_y - ix - 1;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > (ix + 1)) {
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
            for (DefInt ix = 0; ix < region_length; ++ix) {
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
                        if (index_min > (ix + 1)) {
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
                if (index_min > (iy + 1)) {
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
            sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                |domain_min_n_level.at(kYIndex);
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
        } else {
            index_min = 0;
            bool_not_y_max = false;
        }
    }
    if (bool_not_y_max) {
        for (DefInt iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy) * total_length + region_length;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > (ix + 1)) {
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
                for (DefInt ix = 0; ix < region_length; ++ix) {
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
                            if (index_min > (ix + 1)) {
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
                    if (index_min > (iy + 1)) {
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
 * @brief function to find nodes in the region of a given length.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length the length of the region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 * @param[out] ptr_sfbitset_node_overlap pointer to nodes overlapped with periodic boundary.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid.
 */
// for periodic boundaries, nodes along the boundary are assumed to be overlapped where min == max
// example for a length 2 region, where o is the input, * is inner nodes,
// i is minimum periodic boundary, and a is maximum periodic boundary
// i * * * * a
// i * * * * a
// i * * * o a
// i * * * * a
// i * * * * a
// ptr_sfbitset_nodes will include nodes (* and a) and ptr_sfbitset_node_overlap will include nodes (i)
DefInt SFBitsetAux2D::FindNodesInPeriodicRegionCornerOverlap(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes,
    std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, kInvalidSFbitset);
    ptr_sfbitset_node_overlap->clear();
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefSFBitset sfbitset_overlap_y = kInvalidSFbitset, sfbitset_overlap_x = kInvalidSFbitset;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_overlap_y;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy < region_length; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        if (sfbitset_overlap_y != kInvalidSFbitset) {
            sfbitset_overlap_x = sfbitset_overlap_y;
        } else {
            sfbitset_overlap_x = kInvalidSFbitset;
        }
        vec_index_y = (region_length - iy - 1) * total_length + region_length;
        for (DefInt ix = 0; ix < region_length; ++ix) {
            vec_index_x = vec_index_y - ix - 1;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if (sfbitset_overlap_x != kInvalidSFbitset) {
                ptr_sfbitset_node_overlap->emplace_back(
                    std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_x));
            }
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                            vec_index_overlap_y - ix - 1, sfbitset_overlap_x));
                    }
                } else {
                    if (index_min > (ix + 1)) {
                        index_min = ix + 1;
                    }
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            if (sfbitset_overlap_x != kInvalidSFbitset) {
                sfbitset_overlap_x = FindXNeg(sfbitset_overlap_x);
            }
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        if (sfbitset_overlap_y != kInvalidSFbitset) {
            sfbitset_overlap_x = sfbitset_overlap_y;
        } else {
            sfbitset_overlap_x = kInvalidSFbitset;
        }
        vec_index_y = (region_length - iy - 1) * total_length + region_length;
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_x));
                }
            } else {
                index_min = 0;
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefInt ix = 0; ix < region_length; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                vec_index_x = vec_index_y + ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = FindXPos(sfbitset_overlap_x);
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (sfbitset_overlap_x != kInvalidSFbitset) {
                            sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                        }
                    } else {
                        if (index_min > (ix + 1)) {
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
                sfbitset_overlap_y = sfbitset_tmp_y;
                vec_index_overlap_y = vec_index_y;
            } else {
                if (index_min > (iy + 1)) {
                    index_min = iy + 1;
                }
                break;
            }
        } else {
            sfbitset_overlap_y = kInvalidSFbitset;
        }
        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    bool_not_y_max = true;
    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kYIndex)) {
        if (periodic_max.at(kYIndex)) {
            sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                |domain_min_n_level.at(kYIndex);
            sfbitset_overlap_y = sfbitset_tmp_y;
            vec_index_overlap_y = vec_index_y;
        } else {
            index_min = 0;
            bool_not_y_max = false;
        }
    } else {
        sfbitset_overlap_y = kInvalidSFbitset;
    }
    if (bool_not_y_max) {
        for (DefInt iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            if (sfbitset_overlap_y != kInvalidSFbitset) {
                sfbitset_overlap_x = sfbitset_overlap_y;
            } else {
                sfbitset_overlap_x = kInvalidSFbitset;
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy) * total_length + region_length;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_x));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (sfbitset_overlap_x != kInvalidSFbitset) {
                            sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_x));
                        }
                    } else {
                        if (index_min > (ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = FindXNeg(sfbitset_overlap_x);
                }
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            if (sfbitset_overlap_y != kInvalidSFbitset) {
                sfbitset_overlap_x = sfbitset_overlap_y;
            } else {
                sfbitset_overlap_x = kInvalidSFbitset;
            }
            vec_index_y = (region_length + iy) * total_length + region_length;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_x));
                    }
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = FindXPos(sfbitset_overlap_x);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (sfbitset_overlap_x != kInvalidSFbitset) {
                                sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                            }
                        } else {
                            if (index_min > (ix + 1)) {
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
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_y;
                } else {
                    if (index_min > (iy + 1)) {
                        index_min = iy + 1;
                    }
                    break;
                }
            } else {
                sfbitset_overlap_y = kInvalidSFbitset;
            }
        }
    }
    return index_min;
}
DefInt SFBitsetAux2D::FindNodesInPeriodicRegionCenterOverlap(const DefSFBitset& sfbitset_in,
    const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes,
    std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
                + ") or maximum (" + std::to_string(domain_max_n_level.size())
                + ") boundary is not equal to space filling code dimension (2)");
    }
    if (region_length_neg.size() < 2 || region_length_pos.size() < 2) {
        LogManager::LogError("dimension of x (" + std::to_string(region_length_neg.size())
                + ") or y (" + std::to_string(region_length_pos.size())
                + ") search distance is less than space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length_x = region_length_neg[kXIndex] + region_length_pos[kXIndex] + 1;
    DefAmrLUint total_length_y = region_length_neg[kYIndex] + region_length_pos[kYIndex] + 1;
    ptr_sfbitset_nodes->resize(total_length_x * total_length_y);
    ptr_sfbitset_nodes->assign(total_length_x * total_length_y, kInvalidSFbitset);
    ptr_sfbitset_node_overlap->clear();
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefSFBitset sfbitset_overlap_y = kInvalidSFbitset, sfbitset_overlap_x = kInvalidSFbitset;
    DefAmrLUint vec_index_x, vec_index_y,
        vec_index_overlap_y = region_length_neg[kYIndex]* total_length_x + region_length_neg[kXIndex];
    DefInt index_min = region_length_neg[kXIndex] + 1;
    if (index_min> region_length_pos[kXIndex]) {
        index_min = region_length_pos[kXIndex];
    }
    if (index_min> region_length_neg[kYIndex] + 1) {
        index_min = region_length_neg[kYIndex] + 1;
    }
    if (index_min> region_length_pos[kYIndex]) {
        index_min = region_length_pos[kYIndex];
    }
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        if (sfbitset_overlap_y != kInvalidSFbitset) {
            sfbitset_overlap_x = sfbitset_overlap_y;
        } else {
            sfbitset_overlap_x = kInvalidSFbitset;
        }
        vec_index_y = (region_length_neg[kYIndex] - iy) * total_length_x + region_length_neg[kXIndex];
        for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
            vec_index_x = vec_index_y - ix;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if (sfbitset_overlap_x != kInvalidSFbitset) {
                ptr_sfbitset_node_overlap->emplace_back(
                    std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_x));
            }
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                            vec_index_overlap_y - ix, sfbitset_overlap_x));
                    }
                } else {
                    if (index_min > ix) {
                        index_min = ix;
                    }
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            if (sfbitset_overlap_x != kInvalidSFbitset) {
                sfbitset_overlap_x = FindXNeg(sfbitset_overlap_x);
            }
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        if (sfbitset_overlap_y != kInvalidSFbitset) {
            sfbitset_overlap_x = sfbitset_overlap_y;
        } else {
            sfbitset_overlap_x = kInvalidSFbitset;
        }
        vec_index_y = (region_length_neg[kYIndex] - iy) * total_length_x + region_length_neg[kXIndex];
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_x));
                }
            } else {
                index_min = 0;
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                vec_index_x = vec_index_y + ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = FindXPos(sfbitset_overlap_x);
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (sfbitset_overlap_x != kInvalidSFbitset) {
                            sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                        }
                    } else {
                        if (index_min > ix) {
                            index_min = ix;
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
                sfbitset_overlap_y = sfbitset_tmp_y;
                vec_index_overlap_y = vec_index_y;
            } else {
                if (index_min > iy) {
                    index_min = iy;
                }
                break;
            }
        } else {
            sfbitset_overlap_y = kInvalidSFbitset;
        }
        sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
    }
    // positive y direction
    sfbitset_tmp_y = sfbitset_in;
    bool_not_y_max = true;
    if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kYIndex)) {
        if (periodic_max.at(kYIndex)) {
            sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                |domain_min_n_level.at(kYIndex);
            sfbitset_overlap_y = sfbitset_tmp_y;
            vec_index_overlap_y = vec_index_y;
        } else {
            index_min = 0;
            bool_not_y_max = false;
        }
    } else {
        sfbitset_overlap_y = kInvalidSFbitset;
    }
    if (bool_not_y_max) {
        for (DefInt iy = 1; iy < region_length_pos[kYIndex] + 1; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            if (sfbitset_overlap_y != kInvalidSFbitset) {
                sfbitset_overlap_x = sfbitset_overlap_y;
            } else {
                sfbitset_overlap_x = kInvalidSFbitset;
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length_neg[kYIndex] + iy) * total_length_x + region_length_neg[kXIndex];
            for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_x));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (sfbitset_overlap_x != kInvalidSFbitset) {
                            sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_x));
                        }
                    } else {
                        if (index_min > ix) {
                            index_min = ix;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                if (sfbitset_overlap_x != kInvalidSFbitset) {
                    sfbitset_overlap_x = FindXNeg(sfbitset_overlap_x);
                }
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            if (sfbitset_overlap_y != kInvalidSFbitset) {
                sfbitset_overlap_x = sfbitset_overlap_y;
            } else {
                sfbitset_overlap_x = kInvalidSFbitset;
            }
            vec_index_y = (region_length_neg[kYIndex] + iy) * total_length_x + region_length_neg[kXIndex];
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_x));
                    }
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (sfbitset_overlap_x != kInvalidSFbitset) {
                        sfbitset_overlap_x = FindXPos(sfbitset_overlap_x);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (sfbitset_overlap_x != kInvalidSFbitset) {
                                sfbitset_overlap_x = (sfbitset_overlap_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_x));
                            }
                        } else {
                            if (index_min > ix) {
                                index_min = ix;
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
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_y;
                } else {
                    if (index_min > iy) {
                        index_min = iy;
                    }
                    break;
                }
            } else {
                sfbitset_overlap_y = kInvalidSFbitset;
            }
        }
    }
    return index_min;
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length the length of the region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid, otherwise may be undefined.
 */
// example for a length 2 region, o is the input, * is nodes around it
// * * * * *
// * * * * *
// * * o * *
// * * * * *
// * * * * *
DefInt SFBitsetAux2D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy <= region_length; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length - iy) * total_length + region_length;
        for (DefInt ix = 0; ix <= region_length; ++ix) {
            vec_index_x = vec_index_y - ix;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > (ix)) {
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
            for (DefInt ix = 0; ix < region_length; ++ix) {
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
                        if (index_min > (ix + 1)) {
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
                if (index_min > (iy)) {
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
        for (DefInt iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length + iy + 1) * total_length + region_length;
            for (DefInt ix = 0; ix <= region_length; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    if (index_min > (ix)) {
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
                for (DefInt ix = 0; ix < region_length; ++ix) {
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
                            if (index_min > (ix + 1)) {
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
                    if (index_min > (iy + 1)) {
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
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length_neg the length of searching region in negative directions.
 * @param[in] region_length_pos the length of searching region in positive directions.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 */
void SFBitsetAux2D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const DefInt region_length_neg, const DefInt region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = region_length_neg + region_length_pos + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy <= region_length_neg; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length_neg - iy) * total_length + region_length_neg;
        for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
            vec_index_x = vec_index_y - ix;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length_neg - iy) * total_length + region_length_neg;
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            } else {
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
            bool_not_y_max = false;
        }
    }
    if (bool_not_y_max) {
        for (DefInt iy = 0; iy < region_length_pos; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length_neg + iy + 1) * total_length + region_length_neg;
            for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length_neg + iy + 1) * total_length + region_length_neg;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
                    break;
                }
            }
        }
    }
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length_neg lengths of searching region in negative directions.
 * @param[in] region_length_pos lengths of searching region in positive directions.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 */
void SFBitsetAux2D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 2 || domain_max_n_level.size() != 2) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (2)");
    }
    if (region_length_neg.size() < 2 || region_length_pos.size() < 2) {
        LogManager::LogError("size of input searching length should not be less than the dimension");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length_x = region_length_neg[kXIndex] + region_length_pos[kXIndex] + 1;
    DefAmrLUint total_length_y = region_length_neg[kYIndex] + region_length_pos[kYIndex] + 1;
    ptr_sfbitset_nodes->resize(total_length_x * total_length_y);
    ptr_sfbitset_nodes->assign(total_length_x * total_length_y, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_y = sfbitset_in, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y;
    bool bool_not_x_max, bool_not_y_max;
    // negative y direction
    for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length_neg[kYIndex] - iy) * total_length_x + region_length_neg[kXIndex];
        for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
            vec_index_x = vec_index_y - ix;
            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
        }
        sfbitset_tmp_x = sfbitset_tmp_y;
        vec_index_y = (region_length_neg[kYIndex] - iy) * total_length_x + region_length_neg[kXIndex];
        bool_not_x_max = true;
        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
            if (periodic_max.at(kXIndex)) {
                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kXIndex);
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            } else {
                bool_not_x_max = false;
            }
        }
        if (bool_not_x_max) {
            for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
            bool_not_y_max = false;
        }
    }
    if (bool_not_y_max) {
        for (DefInt iy = 0; iy < region_length_pos[kYIndex]; ++iy) {
            sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length_neg[kYIndex] + iy + 1) * total_length_x + region_length_neg[kXIndex];
            for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_min_n_level.at(kXIndex)) {
                if (periodic_min.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                } else {
                    break;
                }
            }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (region_length_neg[kYIndex] + iy + 1) * total_length_x + region_length_neg[kXIndex];
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
                    break;
                }
            }
        }
    }
}
/**
 * @brief function to set index of the input spacing fill code as the given one.
 * @param[in] i_dir the given direction.
 * @param[in] sfbitset_in input spacing fill code at current level.
 * @param[in] sfbitset_target target spacing fill code at current level.
 * @return the new spacing fill code with the index set to the given one.
 */
DefSFBitset SFBitsetAux2D::SetIindexinGivenDirection(const DefInt i_dir, const DefSFBitset& sfbitset_in,
    const DefSFBitset& sfbitset_target) const {
    if (i_dir == kXIndex) {
        return (sfbitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))
            |(sfbitset_target&k0SFBitsetTakeXRef_.at(kRefCurrent_));
    } else if (i_dir == kYIndex) {
        return (sfbitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))
            |(sfbitset_target&k0SFBitsetTakeYRef_.at(kRefCurrent_));
    } else {
        LogManager::LogError("i_dir should be 0 or 1 in SFBitsetAux2D::SetIindexinGivenDirection in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return kInvalidSFbitset;
    }
}
/**
 * @brief function to check if node at current level should exist based on lower level nodes.
 * @param[in] sfbitset_in input node at current level.
 * @param[in] exist_nodes_lower_level existing nodes at one lower level.
 * @return true if the node at current level should exist.
 */
bool SFBitsetAux2D::CheckExistenceCurrentLevel(
    const DefSFBitset& sfbitset_in, const DefMap<DefInt>& exist_nodes_lower_level) const {
    DefSFBitset sfbitset_in_lower = SFBitsetToOneLowerLevel(sfbitset_in);
    if (exist_nodes_lower_level.find(sfbitset_in_lower) == exist_nodes_lower_level.end()) {
        return false;
    }
    DefSFBitset current_bit = sfbitset_in&k0SfBitsetCurrentLevelBits_;
    if (current_bit != 0) {
        if ((current_bit&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != 0) {
            DefSFBitset sfbitset0 = FindXPos(sfbitset_in_lower);
            if (exist_nodes_lower_level.find(sfbitset0)
                == exist_nodes_lower_level.end()) {
                return false;
            }
            if ((current_bit&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != 0) {
                if (exist_nodes_lower_level.find(FindYPos(sfbitset0))
                    == exist_nodes_lower_level.end()) {
                    return false;
                }
            }
        }
        if ((current_bit&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != 0) {
            if (exist_nodes_lower_level.find(FindYPos(sfbitset_in_lower))
                == exist_nodes_lower_level.end()) {
                return false;
            }
        }
    }
    return true;
}
/**
 * @brief function to check if the node is on a given boundary.
 * @param[in] i_dir the given direction.
 * @param[in] sfbitset_in space filling code of a node.
 * @param[in] sfbitset_boundary space filling code of the boundary in the given direction.
 * @return true if the node is on the given boundary.
 */
bool SFBitsetAux2D::CheckIfOnGivenBoundary(const DefInt i_dir,
    const DefSFBitset& sfbitset_in, const DefSFBitset& sfbitset_boundary) const {
    switch (i_dir) {
    case kXIndex:
        if ((sfbitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) == sfbitset_boundary) {
            return true;
        } else {
            return false;
        }
    case kYIndex:
        if ((sfbitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) == sfbitset_boundary) {
            return true;
        } else {
            return false;
        }
    default:
        return false;
    }
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 2
 */
void SFBitsetAux2D::GetMinM1AtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_min,
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
void SFBitsetAux2D::GetMaxP1AtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_max,
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
void SFBitsetAux2D::GetMinAtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_bitsets) const {
    ptr_min_bitsets->resize(2);
    ptr_min_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_min[kXIndex], 0}));
    ptr_min_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_min[kYIndex]}));
}
/**
 * @brief function to calculate spacing fill code of maximum indices at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_max maximum indicies of the computational domain 
 * @param[out] ptr_max_bitsets a pointer to maximum indices
 */
void SFBitsetAux2D::GetMaxAtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_bitsets) const {
    ptr_max_bitsets->resize(2);
    ptr_max_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0}));
    ptr_max_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex]}));
}
/**
* @brief   function to find space filling code of nodes on a cell (2D) which may across two refinement levels
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  map_node_exist nodes at fine grid with i_level - 1 spacing filling code
* @param[in]  map_node_exist_coarse  nodes at coarse grid with i_level - 1 spacing filling code
* @param[out] ptr_sfbitsets nodes on the cell 
* @return  if true indicates the given node belongs to a cell
* @note ptr_sfbitsets[0]:(0, 0); ptr_sfbitsets[1]:(+x, 0);
*       ptr_sfbitsets[2]:(0, +y); ptr_sfbitsets[3]:(+x, +y);
*       the int value of the pair indicating if the nodes is in map_node_exist or map_node_exist_coarse
*/
bool SFBitsetAux2D::SFBitsetBelongToOneCellAcrossTwoLevels(const DefSFBitset& sfbitset_in,
    const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_node_exist_coarse,
    std::array<std::pair<DefSFBitset, DefInt>, 4>* const ptr_sfbitsets) const {
    if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(sfbitset_in) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(sfbitset_in) == 0) {
            return false;
        } else {
            ptr_sfbitsets->at(0).second = 2;
        }
    } else {
        ptr_sfbitsets->at(0).second = 1;
    }
    ptr_sfbitsets->at(0).first = sfbitset_in;
    // (+x, 0)
    ptr_sfbitsets->at(1).first = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(1).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(1).first) == 0) {
            return false;
        } else {
            ptr_sfbitsets->at(1).second = 2;
        }
    } else {
        ptr_sfbitsets->at(1).second = 1;
    }

    // (+x, +y)
    ptr_sfbitsets->at(3).first  = FindYPos(ptr_sfbitsets->at(1).first);
    if (map_node_exist.find(ptr_sfbitsets->at(3).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(3).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(3).first) == 0) {
            return false;
        } else {
            ptr_sfbitsets->at(3).second = 2;
        }
    } else {
        ptr_sfbitsets->at(3).second = 1;
    }
    // (0, +y)
    ptr_sfbitsets->at(2).first = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(2).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(2).first) == 0) {
            return false;
        } else {
            ptr_sfbitsets->at(2).second = 2;
        }
    } else {
        ptr_sfbitsets->at(2).second = 1;
    }
    return true;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
void SFBitsetAux3D::SetSpaceBackground(const std::vector<DefReal>& space_background) {
    if (space_background.size() != 3) {
        LogManager::LogError("dimension of space background should be 3");
    }
    k0SpaceBackground_ = space_background;
}
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
* @param[in]  level_diff  different between target and given levels of nodes.
* @param[in]  sfbitset_in  sfbitset at the lower corner.
* @param[out]  ptr_vec_sfbitsets_higher_level
*              sfbitsets at the higher refinement level in a cell
*              at the lower level.
*/
void SFBitsetAux3D::SFBitsetHigherLevelInACell(
    const DefInt level_diff, const DefSFBitset& sfbitset_corner,
    std::vector<DefSFBitset>* const ptr_vec_sfbitsets_higher_level) const {
    DefSizet num = TwoPowerN(level_diff);
    ptr_vec_sfbitsets_higher_level->resize(num*num*num);
    DefSFBitset sfbitset_x, sfbitset_y, sfbitset_z = SFBitsetToNHigherLevel(level_diff, sfbitset_corner);
    DefSizet yindex, zindex;
    for (DefSizet iz = 0; iz < num; ++iz) {
        zindex = iz * num * num;
        sfbitset_y = sfbitset_z;
        for (DefSizet iy = 0; iy < num; ++iy) {
            yindex = zindex + iy * num;
            sfbitset_x = sfbitset_y;
            for (DefSizet ix = 0; ix < num; ++ix) {
                ptr_vec_sfbitsets_higher_level->at(yindex + ix) = sfbitset_x;
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
* @brief   function to find all neighboring sfbitset (3D).
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
    DefSFBitset sfbitset_tmp0 = sfbitset_center, sfbitset_tmp1 = sfbitset_center,
     sfbitset_tmp2 = sfbitset_center;

    ptr_bitset_neighbors->at(kNodeIndexX0Y0Z0_) = sfbitset_center;
    // node at (-x, 0, 0)
    sfbitset_tmp0 = bool_x_gt0 ? FindXNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Z0_) = sfbitset_tmp0;
    // node at (-x, -y, 0)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZ0_) = sfbitset_tmp1;
    // node at (-x, +y, 0)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZ0_) = sfbitset_tmp1;
    // node at (+x, 0, 0)
    sfbitset_tmp0 = FindXPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Z0_) = sfbitset_tmp0;
    // node at (+x, -y, 0)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZ0_) = sfbitset_tmp1;
    // node at (+x, +y, 0)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZ0_) = sfbitset_tmp1;
    // node at (0, -y, 0)
    sfbitset_tmp0 = bool_y_gt0 ? FindYNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZ0_) = sfbitset_tmp0;
    // node at (0, +y, 0)
    sfbitset_tmp0 = FindYPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZ0_) = sfbitset_tmp0;

    // node at (0, 0, -z)
    sfbitset_tmp0 = bool_z_gt0 ? FindZNeg(sfbitset_center) : sfbitset_center;
    ptr_bitset_neighbors->at(kNodeIndexX0Y0Zn_) = sfbitset_tmp0;
    // node at (-x, 0, -z)
    sfbitset_tmp1 = bool_x_gt0 ? FindXNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Zn_) = sfbitset_tmp1;
    // node at (-x, -y, -z)
    sfbitset_tmp2 = bool_y_gt0 ? FindYNeg(sfbitset_tmp1) : sfbitset_tmp1;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZn_) = sfbitset_tmp2;
    // node at (-x, +y, -z)
    sfbitset_tmp2 = FindYPos(sfbitset_tmp1);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZn_) = sfbitset_tmp2;
    // node at (+x, 0, -z)
    sfbitset_tmp1 = FindXPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Zn_) = sfbitset_tmp1;
    // node at (+x, -y, -z)
    sfbitset_tmp2 = bool_y_gt0 ? FindYNeg(sfbitset_tmp1) : sfbitset_tmp1;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZn_) = sfbitset_tmp2;
    // node at (+x, +y, -z)
    sfbitset_tmp2 = FindYPos(sfbitset_tmp1);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZn_) = sfbitset_tmp2;
    // node at (0, -y, -z)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZn_) = sfbitset_tmp1;
    //  node at (0, +y, -z)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZn_) = sfbitset_tmp1;
    // node at (0, 0, +z)
    sfbitset_tmp0 = FindZPos(sfbitset_center);
    ptr_bitset_neighbors->at(kNodeIndexX0Y0Zp_) = sfbitset_tmp0;
    // node at (-x, 0, +z)
    sfbitset_tmp1 = bool_x_gt0 ? FindXNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexXnY0Zp_) = sfbitset_tmp1;
    // node at (-x, -y, +z)
    sfbitset_tmp2 = bool_y_gt0 ? FindYNeg(sfbitset_tmp1) : sfbitset_tmp1;
    ptr_bitset_neighbors->at(kNodeIndexXnYnZp_) = sfbitset_tmp2;
    // node at (-x, +y, +z)
    sfbitset_tmp2 = FindYPos(sfbitset_tmp1);
    ptr_bitset_neighbors->at(kNodeIndexXnYpZp_) = sfbitset_tmp2;
    // node at (+x, 0, +z)
    sfbitset_tmp1 = FindXPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexXpY0Zp_) = sfbitset_tmp1;
    // node at (+x, -y, +z)
    sfbitset_tmp2 = bool_y_gt0 ? FindYNeg(sfbitset_tmp1) : sfbitset_tmp1;
    ptr_bitset_neighbors->at(kNodeIndexXpYnZp_) = sfbitset_tmp2;
    // node at (+x, +y, +z)
    sfbitset_tmp2 = FindYPos(sfbitset_tmp1);
    ptr_bitset_neighbors->at(kNodeIndexXpYpZp_) = sfbitset_tmp2;
    // node at (0, -y, +z)
    sfbitset_tmp1 = bool_y_gt0 ? FindYNeg(sfbitset_tmp0) : sfbitset_tmp0;
    ptr_bitset_neighbors->at(kNodeIndexX0YnZp_) = sfbitset_tmp1;
    // node at (0, +y, +z)
    sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
    ptr_bitset_neighbors->at(kNodeIndexX0YpZp_) = sfbitset_tmp1;
}
/**
* @brief   function to find all neighboring nodes bonded by lower and upper limits.
* @param[in]  sfbitset_center sfbitset at the center.
* @param[in] bool_periodic_min booleans indicating if the boundary is periodic at maximum domain boundaries.
* @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
* @param[in]  domain_min_n_level   lower bound.
* @param[in]  domain_max_n_level   upper bound.
* @param[out]  ptr_bitset_neighbors   space fill code for neighbors of the give node.
*/
void SFBitsetAux3D::SFBitsetFindAllBondedNeighborsVir(const DefSFBitset& bitset_in,
    const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level, const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_bitset_neighbors) const {
    DefSFBitset sfbitset_tmp0 = bitset_in, sfbitset_tmp1 = bitset_in;
    ptr_bitset_neighbors->clear();
    bool bool_not_x_min = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_min_n_level.at(kXIndex),
        bool_not_x_max = (bitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != domain_max_n_level.at(kXIndex),
        bool_not_y_min = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_min_n_level.at(kYIndex),
        bool_not_y_max = (bitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != domain_max_n_level.at(kYIndex),
        bool_not_z_min = (bitset_in&k0SFBitsetTakeZRef_.at(kRefCurrent_)) != domain_min_n_level.at(kZIndex),
        bool_not_z_max = (bitset_in&k0SFBitsetTakeZRef_.at(kRefCurrent_)) != domain_max_n_level.at(kZIndex);


    if (bool_not_x_min) {
        sfbitset_tmp0 = FindXNeg(bitset_in);
        // node at (-x, 0, 0)
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    } else if (bool_periodic_min.at(kXIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))|domain_max_n_level.at(kXIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_y_min) {
        // node at (-x, -y, 0)
        sfbitset_tmp1 = FindYNeg(sfbitset_tmp0);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    } else if (bool_periodic_min.at(kYIndex)) {
        sfbitset_tmp1 = (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    }
    if (bool_not_z_min) {
        // node at (-x, -y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp1));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (-x, -y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp1));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }
    if (bool_not_y_max) {
        // node at (-x, +y, 0)
        sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    } else if (bool_periodic_max.at(kYIndex)) {
        sfbitset_tmp1 = (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    }
    if (bool_not_z_min) {
        // node at (-x, +y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp1));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (-x, +y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp1));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }
    if (bool_not_z_min) {
        // node at (-x, 0, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (-x, 0, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }

    if (bool_not_x_max) {
        sfbitset_tmp0 = FindXPos(bitset_in);
        // node at (+x, 0, 0)
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    } else if (bool_periodic_max.at(kXIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))|domain_min_n_level.at(kXIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_y_min) {
        // node at (+x, -y, 0)
        sfbitset_tmp1 = FindYNeg(sfbitset_tmp0);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    } else if (bool_periodic_min.at(kYIndex)) {
        sfbitset_tmp1 = (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    }
    if (bool_not_z_min) {
        // node at (+x, -y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp1));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (+x, -y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp1));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }
    if (bool_not_y_max) {
        // node at (+x, +y, 0)
        sfbitset_tmp1 = FindYPos(sfbitset_tmp0);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    } else if (bool_periodic_max.at(kYIndex)) {
        sfbitset_tmp1 = (sfbitset_tmp0&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp1);
    }
    if (bool_not_z_min) {
        // node at (+x, +y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp1));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (+x, +y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp1));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp1&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }
    if (bool_not_z_min) {
        // node at (+x, 0, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (+x, 0, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }

    if (bool_not_y_min) {
        // node at (0, -y, 0)
        sfbitset_tmp0 = FindYNeg(bitset_in);
        ptr_bitset_neighbors->emplace_back(bitset_in);
    } else if (bool_periodic_min.at(kYIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_max_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_z_min) {
        // node at (0, -y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (0, -y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }

    if (bool_not_y_max) {
        // node at (0, +y, 0)
        sfbitset_tmp0 = FindYPos(bitset_in);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    } else if (bool_periodic_max.at(kYIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))|domain_min_n_level.at(kYIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_z_min) {
        // node at (0, +y, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(sfbitset_tmp0));
    } else if (bool_periodic_min.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex));
    }
    if (bool_not_z_max) {
        // node at (0, +y, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(sfbitset_tmp0));
    } else if (bool_periodic_max.at(kZIndex)) {
        ptr_bitset_neighbors->emplace_back(
            (sfbitset_tmp0&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex));
    }

    if (bool_not_z_min) {
        // node at (0, 0, -z)
        ptr_bitset_neighbors->emplace_back(FindZNeg(bitset_in));
    } else if (bool_periodic_min.at(kZIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_max_n_level.at(kZIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
    }
    if (bool_not_z_max) {
        // node at (0, 0, +z)
        ptr_bitset_neighbors->emplace_back(FindZPos(bitset_in));
    } else if (bool_periodic_max.at(kZIndex)) {
        sfbitset_tmp0 = (bitset_in&k0SFBitsetTakeZRef_.at(kRefOthers_))|domain_min_n_level.at(kZIndex);
        ptr_bitset_neighbors->emplace_back(sfbitset_tmp0);
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
    const DefInt region_length, const std::vector<DefSFBitset>& domain_min_m1_n_level,
    const std::vector<DefSFBitset>& domain_max_p1_n_level, std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_m1_n_level.size() != 3 || domain_max_p1_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_m1_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_p1_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID
    DefSFBitset sfbitset_tmp_y, sfbitset_tmp_x, sfbitset_tmp_z = sfbitset_in;
    ptr_sfbitset_nodes->clear();
    // negative z direction
    for (DefInt iz = 0; iz <= region_length; ++iz) {
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
         != domain_min_m1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefInt iy = 0; iy <= region_length; ++iy) {
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
    for (DefInt iz = 0; iz < region_length; ++iz) {
        sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
         != domain_max_p1_n_level.at(kZIndex)) {
            sfbitset_tmp_y = sfbitset_tmp_z;
            // negative y direction
            for (DefInt iy = 0; iy <= region_length; ++iy) {
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_min_m1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                 != domain_max_p1_n_level.at(kYIndex)) {
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix <= region_length; ++ix) {
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_min_m1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                         != domain_max_p1_n_level.at(kXIndex)) {
                            ptr_sfbitset_nodes->emplace_back(sfbitset_tmp_x);
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
DefInt SFBitsetAux3D::FindNodesInPeriodicRegionCorner(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefInt iz = 0; iz < region_length; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length - iz - 1) * total_length  + region_length;
        for (DefInt iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > (ix + 1)) {
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
                for (DefInt ix = 0; ix < region_length; ++ix) {
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
                            if (index_min > (ix + 1)) {
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
                    if (index_min > (iy + 1)) {
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
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy) * total_length + region_length;
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > (ix + 1)) {
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
                    for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                if (index_min > (ix + 1)) {
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
                        if (index_min > (iy + 1)) {
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
                if (index_min > (iz + 1)) {
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
            sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                |domain_min_n_level.at(kZIndex);
            sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
        } else {
            index_min = 0;
            bool_not_z_max = false;
        }
    }
    if (bool_not_z_max) {
        for (DefInt iz = 0; iz < region_length; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length + iz) * total_length + region_length;
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > (ix + 1)) {
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
                    for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                if (index_min > (ix + 1)) {
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
                        if (index_min > (iy + 1)) {
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
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefInt iy = 0; iy < region_length; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy) * total_length + region_length;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        vec_index_x = vec_index_y - ix - 1;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > (ix + 1)) {
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
                        for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                    if (index_min > (ix + 1)) {
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
                            if (index_min > (iy + 1)) {
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
                    if (index_min > (iz + 1)) {
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
 * @brief function to find nodes in the region of a given length.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length the length of the region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 * @param[out] ptr_sfbitset_node_overlap pointer to nodes overlapped with periodic boundary.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid.
 */
DefInt SFBitsetAux3D::FindNodesInPeriodicRegionCornerOverlap(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes,
    std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, kInvalidSFbitset);
    ptr_sfbitset_node_overlap->clear();
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    DefSFBitset sfbitset_overlap_yx = kInvalidSFbitset, sfbitset_overlap_zx_y = kInvalidSFbitset,
        sfbitset_overlap_zx = kInvalidSFbitset,
        sfbitset_overlap_zyx = kInvalidSFbitset, sfbitset_overlap_z = kInvalidSFbitset,
        sfbitset_overlap_y = kInvalidSFbitset, sfbitset_overlap_zy = kInvalidSFbitset;
    DefAmrLUint vec_index_overlap_z, vec_index_overlap_zny, vec_index_overlap_zy, vec_index_overlap_y;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    bool periodic_z = false, periodic_y = false;
    // negative z direction
    for (DefInt iz = 0; iz < region_length; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length - iz - 1) * total_length  + region_length;
        if (periodic_z) {
            sfbitset_overlap_zy = sfbitset_overlap_z;
            sfbitset_overlap_zx_y = sfbitset_overlap_z;
        } else {
            sfbitset_overlap_zy = kInvalidSFbitset;
            sfbitset_overlap_zx_y = kInvalidSFbitset;
        }
        for (DefInt iy = 0; iy < region_length; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
            if (periodic_z) {
                vec_index_overlap_zy = (vec_index_overlap_z - iy - 1) * total_length + region_length;
            }
            sfbitset_overlap_yx = sfbitset_overlap_y;
            sfbitset_overlap_zx = sfbitset_overlap_zx_y;
            sfbitset_overlap_zyx = sfbitset_overlap_zy;
            for (DefInt ix = 0; ix < region_length; ++ix) {
                vec_index_x = vec_index_y - ix - 1;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (periodic_z) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                }
                if (periodic_y) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                }
                if (periodic_z && periodic_y) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                        }
                    } else {
                        if (index_min > (ix + 1)) {
                            index_min = ix + 1;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                if (periodic_z) {
                    sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                }
                if (periodic_y) {
                    sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                }
                if (periodic_z&&periodic_y) {
                    sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                }
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            sfbitset_overlap_yx = sfbitset_overlap_y;
            sfbitset_overlap_zx = sfbitset_overlap_zx_y;
            sfbitset_overlap_zyx = sfbitset_overlap_zy;
            vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
            if (periodic_z) {
                vec_index_overlap_zy = (vec_index_overlap_z - iy - 1) * total_length + region_length;
            }
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                    if (periodic_z) {
                        sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy + region_length - 1, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny + region_length - 1, sfbitset_overlap_zyx));
                    }
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > (ix + 1)) {
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
                    sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    periodic_y = true;
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_y;
                    vec_index_overlap_zny = (vec_index_overlap_z - iy - 1) * total_length + region_length;
                } else {
                    if (index_min > (iy + 1)) {
                        index_min = iy + 1;
                    }
                    break;
                }
            } else {
                periodic_y = false;
                sfbitset_overlap_y = kInvalidSFbitset;
                sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                    &kInvalidSFbitset);
            }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            if (periodic_z) {
                sfbitset_overlap_zx_y = FindYNeg(sfbitset_overlap_zx_y);
            }
        }
        // positive y direction
        sfbitset_tmp_y = sfbitset_tmp_z;
        if (periodic_z) {
            sfbitset_overlap_zy = sfbitset_overlap_z;
            sfbitset_overlap_zx_y = sfbitset_overlap_z;
        } else {
            sfbitset_overlap_zy = kInvalidSFbitset;
            sfbitset_overlap_zx_y = kInvalidSFbitset;
        }
        bool_not_y_max = true;
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kYIndex)) {
            if (periodic_max.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                periodic_y = true;
                sfbitset_overlap_y = sfbitset_tmp_y;
                vec_index_overlap_y = vec_index_z * total_length + region_length;
                vec_index_overlap_zny = (vec_index_overlap_z) * total_length + region_length;
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        } else {
            periodic_y = false;
            sfbitset_overlap_y = kInvalidSFbitset;
            sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                &kInvalidSFbitset);
        }
        if (bool_not_y_max) {
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if (periodic_z) {
                    sfbitset_overlap_zx_y = FindYPos(sfbitset_overlap_zx_y);
                }
                vec_index_y = (vec_index_z + iy) * total_length + region_length;
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length + region_length;
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > (ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                    }
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                vec_index_y = (vec_index_z + iy) * total_length + region_length;
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length + region_length;
                }
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + region_length - 1, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + region_length - 1, sfbitset_overlap_zyx));
                        }
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > (ix + 1)) {
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
                        sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        periodic_y = true;
                        sfbitset_overlap_y = sfbitset_tmp_y;
                        vec_index_overlap_y = vec_index_y;
                        vec_index_overlap_zny =  (vec_index_overlap_z + iy) * total_length + region_length;
                    } else {
                        if (index_min > (iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                }  else {
                    periodic_y = false;
                    sfbitset_overlap_y = kInvalidSFbitset;
                    sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                        &kInvalidSFbitset);
                }
            }
        }
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
            == domain_min_n_level.at(kZIndex)) {
            if (periodic_min.at(kZIndex)) {
                sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kZIndex);
                periodic_z = true;
                sfbitset_overlap_z = sfbitset_tmp_z;
                vec_index_overlap_z = vec_index_z;
            } else {
                if (index_min > (iz + 1)) {
                    index_min = iz + 1;
                }
                break;
            }
        } else {
            periodic_z = false;
            sfbitset_overlap_z = kInvalidSFbitset;
        }
        sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    bool_not_z_max = true;
    if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kZIndex)) {
        if (periodic_max.at(kZIndex)) {
            sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                |domain_min_n_level.at(kZIndex);
            periodic_z = true;
            sfbitset_overlap_z = sfbitset_tmp_z;
            vec_index_overlap_z = vec_index_z;
        } else {
            index_min = 0;
            bool_not_z_max = false;
        }
    } else {
        periodic_z = false;
        sfbitset_overlap_z = kInvalidSFbitset;
    }

    if (bool_not_z_max) {
        for (DefInt iz = 0; iz < region_length; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            if (periodic_z) {
                sfbitset_overlap_zy = sfbitset_overlap_z;
                sfbitset_overlap_zx_y = sfbitset_overlap_z;
            } else {
                sfbitset_overlap_zy = kInvalidSFbitset;
                sfbitset_overlap_zx_y = kInvalidSFbitset;
            }
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length + iz) * total_length + region_length;
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z - iy - 1) * total_length + region_length;
                }
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                for (DefInt ix = 0; ix < region_length; ++ix) {
                    vec_index_x = vec_index_y - ix - 1;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                    }
                    if (periodic_z && periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > (ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                    }
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                vec_index_y = (vec_index_z - iy - 1) * total_length + region_length;
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z - iy - 1) * total_length + region_length;
                }
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + region_length - 1, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + region_length - 1, sfbitset_overlap_zyx));
                        }
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > (ix + 1)) {
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
                        sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        periodic_y = true;
                        sfbitset_overlap_y = sfbitset_tmp_y;
                        vec_index_overlap_y = vec_index_y;
                        vec_index_overlap_zny = (vec_index_overlap_z - iy - 1) * total_length + region_length;
                    } else {
                        if (index_min > (iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                } else {
                    periodic_y = false;
                    sfbitset_overlap_y = kInvalidSFbitset;
                    sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                        &kInvalidSFbitset);
                }
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                if (periodic_z) {
                    sfbitset_overlap_zx_y = FindYNeg(sfbitset_overlap_zx_y);
                }
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            if (periodic_z) {
                sfbitset_overlap_zy = sfbitset_overlap_z;
                sfbitset_overlap_zx_y = sfbitset_overlap_z;
            } else {
                sfbitset_overlap_zy = kInvalidSFbitset;
                sfbitset_overlap_zx_y = kInvalidSFbitset;
            }
            bool_not_y_max = true;
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    periodic_y = true;
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_z * total_length + region_length;
                    vec_index_overlap_zny = (vec_index_overlap_z) * total_length + region_length;
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            } else {
                periodic_y = false;
                sfbitset_overlap_y = kInvalidSFbitset;
                sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                    &kInvalidSFbitset);
            }
            if (bool_not_y_max) {
                for (DefInt iy = 0; iy < region_length; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    if (periodic_z) {
                        sfbitset_overlap_zx_y = FindYPos(sfbitset_overlap_zx_y);
                    }
                    vec_index_y = (vec_index_z + iy) * total_length + region_length;
                    if (periodic_z) {
                        vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length + region_length;
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    sfbitset_overlap_yx = sfbitset_overlap_y;
                    sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                    sfbitset_overlap_zyx = sfbitset_overlap_zy;
                    for (DefInt ix = 0; ix < region_length; ++ix) {
                        vec_index_x = vec_index_y - ix - 1;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_min_n_level.at(kXIndex)) {
                            if (periodic_min.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy - ix - 1, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y - ix - 1, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny - ix - 1, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > (ix + 1)) {
                                    index_min = ix + 1;
                                }
                                break;
                            }
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                        }
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    sfbitset_overlap_yx = sfbitset_overlap_y;
                    sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                    sfbitset_overlap_zyx = sfbitset_overlap_zy;
                    vec_index_y = (vec_index_z + iy) * total_length + region_length;
                    if (periodic_z) {
                        vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length + region_length;
                    }
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy + region_length - 1, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + region_length - 1, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny + region_length - 1, sfbitset_overlap_zyx));
                            }
                        } else {
                            index_min = 0;
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefInt ix = 0; ix < region_length; ++ix) {
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                            vec_index_x = vec_index_y + ix;
                            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                            if (periodic_z) {
                                sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                            }
                            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                                == domain_max_n_level.at(kXIndex)) {
                                if (periodic_max.at(kXIndex)) {
                                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_x, sfbitset_tmp_x));
                                    if (periodic_z) {
                                        sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                    }
                                    if (periodic_y) {
                                        sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                    }
                                    if (periodic_z&&periodic_y) {
                                        sfbitset_overlap_zyx = (sfbitset_overlap_zyx
                                            &k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                    }
                                } else {
                                    if (index_min > (ix + 1)) {
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
                            sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            periodic_y = true;
                            sfbitset_overlap_y = sfbitset_tmp_y;
                            vec_index_overlap_y = vec_index_y;
                            vec_index_overlap_zny =  (vec_index_overlap_z + iy) * total_length + region_length;
                        } else {
                            if (index_min > (iy + 1)) {
                                index_min = iy + 1;
                            }
                            break;
                        }
                    }  else {
                        periodic_y = false;
                        sfbitset_overlap_y = kInvalidSFbitset;
                        sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                            &kInvalidSFbitset);
                    }
                }
            }
            if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kZIndex)) {
                if (periodic_max.at(kZIndex)) {
                    sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kZIndex);
                    periodic_z = true;
                    sfbitset_overlap_z = sfbitset_tmp_z;
                    vec_index_overlap_z = vec_index_z;
                } else {
                    if (index_min > (iz + 1)) {
                        index_min = iz + 1;
                    }
                    break;
                }
            } else {
                periodic_z = false;
                sfbitset_overlap_z = kInvalidSFbitset;
            }
        }
    }
    return index_min;
}
DefInt SFBitsetAux3D::FindNodesInPeriodicRegionCenterOverlap(const DefSFBitset& sfbitset_in,
    const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes,
    std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
            + ") or maximum (" + std::to_string(domain_max_n_level.size())
            + ") boundary is not equal to space filling code dimension (3)");
    }
    if (region_length_neg.size() < 3 || region_length_pos.size() < 3) {
        LogManager::LogError("dimension of x (" + std::to_string(region_length_neg.size())
            + ") or y (" + std::to_string(region_length_pos.size())
            + ") search distance is less than space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID

    DefAmrLUint total_length_x = region_length_neg[kXIndex] + region_length_pos[kXIndex] + 1;
    DefAmrLUint total_length_y = region_length_neg[kYIndex] + region_length_pos[kYIndex] + 1;
    DefAmrLUint total_length_z = region_length_neg[kZIndex] + region_length_pos[kZIndex] + 1;
    ptr_sfbitset_nodes->resize(total_length_x * total_length_y * total_length_z);
    ptr_sfbitset_nodes->assign(total_length_x * total_length_y * total_length_z, kInvalidSFbitset);
    ptr_sfbitset_node_overlap->clear();
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    DefSFBitset sfbitset_overlap_yx = kInvalidSFbitset, sfbitset_overlap_zx_y = kInvalidSFbitset,
        sfbitset_overlap_zx = kInvalidSFbitset,
        sfbitset_overlap_zyx = kInvalidSFbitset, sfbitset_overlap_z = kInvalidSFbitset,
        sfbitset_overlap_y = kInvalidSFbitset, sfbitset_overlap_zy = kInvalidSFbitset;
    DefAmrLUint vec_index_overlap_z, vec_index_overlap_zny, vec_index_overlap_zy, vec_index_overlap_y;
    DefInt index_min = region_length_neg[kXIndex] + 1;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    bool periodic_z = false, periodic_y = false;
    if (index_min> region_length_pos[kXIndex]) {
        index_min = region_length_pos[kXIndex];
    }
    if (index_min> region_length_neg[kYIndex] + 1) {
        index_min = region_length_neg[kYIndex] + 1;
    }
    if (index_min> region_length_pos[kYIndex]) {
        index_min = region_length_pos[kYIndex];
    }
    if (index_min> region_length_neg[kZIndex] + 1) {
        index_min = region_length_neg[kZIndex] + 1;
    }
    if (index_min> region_length_pos[kZIndex]) {
        index_min = region_length_pos[kZIndex];
    }
    // negative z direction
    for (DefInt iz = 0; iz <= region_length_neg[kZIndex]; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length_neg[kZIndex] - iz) * total_length_y  + region_length_neg[kYIndex];
        if (periodic_z) {
            sfbitset_overlap_zy = sfbitset_overlap_z;
            sfbitset_overlap_zx_y = sfbitset_overlap_z;
        } else {
            sfbitset_overlap_zy = kInvalidSFbitset;
            sfbitset_overlap_zx_y = kInvalidSFbitset;
        }
        for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
            if (periodic_z) {
                vec_index_overlap_zy = (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
            }
            sfbitset_overlap_yx = sfbitset_overlap_y;
            sfbitset_overlap_zx = sfbitset_overlap_zx_y;
            sfbitset_overlap_zyx = sfbitset_overlap_zy;
            for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if (periodic_z) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                }
                if (periodic_y) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                }
                if (periodic_z && periodic_y) {
                    ptr_sfbitset_node_overlap->emplace_back(
                        std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                }
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                        }
                    } else {
                        if (index_min > ix) {
                            index_min = ix;
                        }
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                if (periodic_z) {
                    sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                }
                if (periodic_y) {
                    sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                }
                if (periodic_z&&periodic_y) {
                    sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                }
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            sfbitset_overlap_yx = sfbitset_overlap_y;
            sfbitset_overlap_zx = sfbitset_overlap_zx_y;
            sfbitset_overlap_zyx = sfbitset_overlap_zy;
            vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
            if (periodic_z) {
                vec_index_overlap_zy = (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
            }
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                    if (periodic_z) {
                        sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy + region_length_neg[kXIndex], sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny +  region_length_neg[kXIndex], sfbitset_overlap_zyx));
                    }
                } else {
                    index_min = 0;
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                    sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    vec_index_x = vec_index_y + ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > ix) {
                                index_min = ix;
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
                    sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_max_n_level.at(kYIndex);
                    periodic_y = true;
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_y;
                    vec_index_overlap_zny = (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
                } else {
                    if (index_min > iy) {
                        index_min = iy;
                    }
                    break;
                }
            } else {
                periodic_y = false;
                sfbitset_overlap_y = kInvalidSFbitset;
                sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                    &kInvalidSFbitset);
            }
            sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            if (periodic_z) {
                sfbitset_overlap_zx_y = FindYNeg(sfbitset_overlap_zx_y);
            }
        }
        // positive y direction
        sfbitset_tmp_y = sfbitset_tmp_z;
        if (periodic_z) {
            sfbitset_overlap_zy = sfbitset_overlap_z;
            sfbitset_overlap_zx_y = sfbitset_overlap_z;
        } else {
            sfbitset_overlap_zy = kInvalidSFbitset;
            sfbitset_overlap_zx_y = kInvalidSFbitset;
        }
        bool_not_y_max = true;
        if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
            == domain_max_n_level.at(kYIndex)) {
            if (periodic_max.at(kYIndex)) {
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                periodic_y = true;
                sfbitset_overlap_y = sfbitset_tmp_y;
                vec_index_overlap_y = vec_index_z * total_length_x + region_length_neg[kXIndex];
                vec_index_overlap_zny = vec_index_overlap_z * total_length_x + region_length_neg[kXIndex];
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        } else {
            periodic_y = false;
            sfbitset_overlap_y = kInvalidSFbitset;
            sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                &kInvalidSFbitset);
        }
        if (bool_not_y_max) {
            for (DefInt iy = 1; iy < region_length_pos[kYIndex] + 1; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                if (periodic_z) {
                    sfbitset_overlap_zx_y = FindYPos(sfbitset_overlap_zx_y);
                }
                vec_index_y = (vec_index_z + iy) * total_length_x + region_length_neg[kXIndex];
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                    }
                    if (periodic_z&&periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > ix) {
                                index_min = ix;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                    }
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                vec_index_y = (vec_index_z + iy) * total_length_x + region_length_neg[kXIndex];
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                }
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + region_length_neg[kXIndex], sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(
                                    vec_index_overlap_zny + region_length_neg[kXIndex], sfbitset_overlap_zyx));
                        }
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > ix) {
                                    index_min = ix;
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
                        sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kYIndex);
                        periodic_y = true;
                        sfbitset_overlap_y = sfbitset_tmp_y;
                        vec_index_overlap_y = vec_index_y;
                        vec_index_overlap_zny =
                            (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                    } else {
                        if (index_min > iy) {
                            index_min = iy;
                        }
                        break;
                    }
                }  else {
                    periodic_y = false;
                    sfbitset_overlap_y = kInvalidSFbitset;
                    sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                        &kInvalidSFbitset);
                }
            }
        }
        if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
            == domain_min_n_level.at(kZIndex)) {
            if (periodic_min.at(kZIndex)) {
                sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                    |domain_max_n_level.at(kZIndex);
                periodic_z = true;
                sfbitset_overlap_z = sfbitset_tmp_z;
                vec_index_overlap_z = vec_index_z;
            } else {
                if (index_min > (iz + 1)) {
                    index_min = iz + 1;
                }
                break;
            }
        } else {
            periodic_z = false;
            sfbitset_overlap_z = kInvalidSFbitset;
        }
        sfbitset_tmp_z = FindZNeg(sfbitset_tmp_z);
    }
    // positive z direction
    sfbitset_tmp_z = sfbitset_in;
    bool_not_z_max = true;
    if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
        == domain_max_n_level.at(kZIndex)) {
        if (periodic_max.at(kZIndex)) {
            sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                |domain_min_n_level.at(kZIndex);
            periodic_z = true;
            sfbitset_overlap_z = sfbitset_tmp_z;
            vec_index_overlap_z = vec_index_z;
        } else {
            index_min = 0;
            bool_not_z_max = false;
        }
    } else {
        periodic_z = false;
        sfbitset_overlap_z = kInvalidSFbitset;
    }

    if (bool_not_z_max) {
        for (DefInt iz = 1; iz < region_length_pos[kZIndex] + 1; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            if (periodic_z) {
                sfbitset_overlap_zy = sfbitset_overlap_z;
                sfbitset_overlap_zx_y = sfbitset_overlap_z;
            } else {
                sfbitset_overlap_zy = kInvalidSFbitset;
                sfbitset_overlap_zx_y = kInvalidSFbitset;
            }
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length_neg[kZIndex] + iz) * total_length_y  + region_length_neg[kYIndex];
            for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
                }
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if (periodic_z) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                    }
                    if (periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                    }
                    if (periodic_z && periodic_y) {
                        ptr_sfbitset_node_overlap->emplace_back(
                            std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                    }
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                            }
                        } else {
                            if (index_min > (ix + 1)) {
                                index_min = ix + 1;
                            }
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    if (periodic_z) {
                        sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                    }
                    if (periodic_y) {
                        sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                    }
                    if (periodic_z&&periodic_y) {
                        sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                    }
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                sfbitset_overlap_yx = sfbitset_overlap_y;
                sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                sfbitset_overlap_zyx = sfbitset_overlap_zy;
                vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
                if (periodic_z) {
                    vec_index_overlap_zy = (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
                }
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                        if (periodic_z) {
                            sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + region_length_neg[kXIndex], sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                                vec_index_overlap_zny + region_length_neg[kXIndex], sfbitset_overlap_zyx));
                        }
                    } else {
                        index_min = 0;
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        vec_index_x = vec_index_y + ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_max_n_level.at(kXIndex)) {
                            if (periodic_max.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > ix) {
                                    index_min = ix;
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
                        sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kYIndex);
                        periodic_y = true;
                        sfbitset_overlap_y = sfbitset_tmp_y;
                        vec_index_overlap_y = vec_index_y;
                        vec_index_overlap_zny =
                            (vec_index_overlap_z - iy) * total_length_x + region_length_neg[kXIndex];
                    } else {
                        if (index_min > (iy + 1)) {
                            index_min = iy + 1;
                        }
                        break;
                    }
                } else {
                    periodic_y = false;
                    sfbitset_overlap_y = kInvalidSFbitset;
                    sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                        &kInvalidSFbitset);
                }
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                if (periodic_z) {
                    sfbitset_overlap_zx_y = FindYNeg(sfbitset_overlap_zx_y);
                }
            }
            // positive y direction
            sfbitset_tmp_y = sfbitset_tmp_z;
            if (periodic_z) {
                sfbitset_overlap_zy = sfbitset_overlap_z;
                sfbitset_overlap_zx_y = sfbitset_overlap_z;
            } else {
                sfbitset_overlap_zy = kInvalidSFbitset;
                sfbitset_overlap_zx_y = kInvalidSFbitset;
            }
            bool_not_y_max = true;
            if ((sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kYIndex)) {
                if (periodic_max.at(kYIndex)) {
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    periodic_y = true;
                    sfbitset_overlap_y = sfbitset_tmp_y;
                    vec_index_overlap_y = vec_index_z * total_length_x + region_length_neg[kXIndex];
                    vec_index_overlap_zny = vec_index_overlap_z * total_length_x + region_length_neg[kXIndex];
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            } else {
                periodic_y = false;
                sfbitset_overlap_y = kInvalidSFbitset;
                sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                    &kInvalidSFbitset);
            }
            if (bool_not_y_max) {
                for (DefInt iy = 1; iy < region_length_pos[kYIndex] + 1; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    if (periodic_z) {
                        sfbitset_overlap_zx_y = FindYPos(sfbitset_overlap_zx_y);
                    }
                    vec_index_y = (vec_index_z + iy) * total_length_x + region_length_neg[kXIndex];
                    if (periodic_z) {
                        vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    sfbitset_overlap_yx = sfbitset_overlap_y;
                    sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                    sfbitset_overlap_zyx = sfbitset_overlap_zy;
                    for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                        vec_index_x = vec_index_y - ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if (periodic_z) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                        }
                        if (periodic_y) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                        }
                        if (periodic_z&&periodic_y) {
                            ptr_sfbitset_node_overlap->emplace_back(
                                std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                        }
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                            == domain_min_n_level.at(kXIndex)) {
                            if (periodic_min.at(kXIndex)) {
                                sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_max_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_x, sfbitset_tmp_x));
                                if (periodic_z) {
                                    sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zy - ix, sfbitset_overlap_zx));
                                }
                                if (periodic_y) {
                                    sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_y - ix, sfbitset_overlap_yx));
                                }
                                if (periodic_z&&periodic_y) {
                                    sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_max_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_overlap_zny - ix, sfbitset_overlap_zyx));
                                }
                            } else {
                                if (index_min > ix) {
                                    index_min = ix;
                                }
                                break;
                            }
                        }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        if (periodic_z) {
                            sfbitset_overlap_zx = FindXNeg(sfbitset_overlap_zx);
                        }
                        if (periodic_y) {
                            sfbitset_overlap_yx = FindXNeg(sfbitset_overlap_yx);
                        }
                        if (periodic_z&&periodic_y) {
                            sfbitset_overlap_zyx = FindXNeg(sfbitset_overlap_zyx);
                        }
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    sfbitset_overlap_yx = sfbitset_overlap_y;
                    sfbitset_overlap_zx = sfbitset_overlap_zx_y;
                    sfbitset_overlap_zyx = sfbitset_overlap_zy;
                    vec_index_y = (vec_index_z + iy) * total_length_x + region_length_neg[kXIndex];
                    if (periodic_z) {
                        vec_index_overlap_zy = (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                    }
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            ptr_sfbitset_node_overlap->emplace_back(std::make_pair(vec_index_y, sfbitset_tmp_x));
                            if (periodic_z) {
                                sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                                    vec_index_overlap_zy + region_length_neg[kXIndex], sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                                    vec_index_overlap_y + region_length_neg[kXIndex], sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = (sfbitset_overlap_zyx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                    |domain_min_n_level.at(kXIndex);
                                ptr_sfbitset_node_overlap->emplace_back(std::make_pair(
                                    vec_index_overlap_zny + region_length_neg[kXIndex], sfbitset_overlap_zyx));
                            }
                        } else {
                            index_min = 0;
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefInt ix = 1; ix < region_length_pos[kXIndex] + 1; ++ix) {
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                            vec_index_x = vec_index_y + ix;
                            ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                            if (periodic_z) {
                                sfbitset_overlap_zx = FindXPos(sfbitset_overlap_zx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                            }
                            if (periodic_y) {
                                sfbitset_overlap_yx = FindXPos(sfbitset_overlap_yx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                            }
                            if (periodic_z&&periodic_y) {
                                sfbitset_overlap_zyx = FindXPos(sfbitset_overlap_zyx);
                                ptr_sfbitset_node_overlap->emplace_back(
                                    std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                            }
                            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                                == domain_max_n_level.at(kXIndex)) {
                                if (periodic_max.at(kXIndex)) {
                                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                        |domain_min_n_level.at(kXIndex);
                                    ptr_sfbitset_node_overlap->emplace_back(
                                        std::make_pair(vec_index_x, sfbitset_tmp_x));
                                    if (periodic_z) {
                                        sfbitset_overlap_zx = (sfbitset_overlap_zx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_zy + ix, sfbitset_overlap_zx));
                                    }
                                    if (periodic_y) {
                                        sfbitset_overlap_yx = (sfbitset_overlap_yx&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_y + ix, sfbitset_overlap_yx));
                                    }
                                    if (periodic_z&&periodic_y) {
                                        sfbitset_overlap_zyx = (sfbitset_overlap_zyx
                                            &k0SFBitsetTakeXRef_.at(kRefOthers_))
                                            |domain_min_n_level.at(kXIndex);
                                        ptr_sfbitset_node_overlap->emplace_back(
                                            std::make_pair(vec_index_overlap_zny + ix, sfbitset_overlap_zyx));
                                    }
                                } else {
                                    if (index_min > ix) {
                                        index_min = ix;
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
                            sfbitset_overlap_zy = (sfbitset_overlap_zy&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            sfbitset_overlap_zx_y = (sfbitset_overlap_zx_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kYIndex);
                            periodic_y = true;
                            sfbitset_overlap_y = sfbitset_tmp_y;
                            vec_index_overlap_y = vec_index_y;
                            vec_index_overlap_zny =
                                (vec_index_overlap_z + iy) * total_length_x + region_length_neg[kXIndex];
                        } else {
                            if (index_min > iy) {
                                index_min = iy;
                            }
                            break;
                        }
                    }  else {
                        periodic_y = false;
                        sfbitset_overlap_y = kInvalidSFbitset;
                        sfbitset_overlap_zy = sfbitset_overlap_zy|(k0SFBitsetTakeYRef_.at(kRefCurrent_)
                            &kInvalidSFbitset);
                    }
                }
            }
            if ((sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kZIndex)) {
                if (periodic_max.at(kZIndex)) {
                    sfbitset_tmp_z = (sfbitset_tmp_z&k0SFBitsetTakeZRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kZIndex);
                    periodic_z = true;
                    sfbitset_overlap_z = sfbitset_tmp_z;
                    vec_index_overlap_z = vec_index_z;
                } else {
                    if (index_min > iz) {
                        index_min = iz;
                    }
                    break;
                }
            } else {
                periodic_z = false;
                sfbitset_overlap_z = kInvalidSFbitset;
            }
        }
    }
    return index_min;
}
/**
 * @brief function to set index of the input spacing fill code as the given one.
 * @param[in] i_dir the given direction.
 * @param[in] sfbitset_in input spacing fill code at current level.
 * @param[in] sfbitset_target target spacing fill code at current level.
 * @return the new spacing fill code with the index set to the given one.
 */
DefSFBitset SFBitsetAux3D::SetIindexinGivenDirection(const DefInt i_dir, const DefSFBitset& sfbitset_in,
    const DefSFBitset& sfbitset_target) const {
    if (i_dir == kXIndex) {
        return (sfbitset_in&k0SFBitsetTakeXRef_.at(kRefOthers_))
            |(sfbitset_target&k0SFBitsetTakeXRef_.at(kRefCurrent_));
    } else if (i_dir == kYIndex) {
        return (sfbitset_in&k0SFBitsetTakeYRef_.at(kRefOthers_))
            |(sfbitset_target&k0SFBitsetTakeYRef_.at(kRefCurrent_));
    } else if (i_dir == kZIndex) {
        return (sfbitset_in&k0SFBitsetTakeZRef_.at(kRefOthers_))
            |(sfbitset_target&k0SFBitsetTakeZRef_.at(kRefCurrent_));
    } else {
        LogManager::LogError("invalid direction index " + std::to_string(i_dir)
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return kInvalidSFbitset;
    }
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length the length of the region.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 * @return number of node in each direction within the region and domain range.
 * @note only indics within the return values are valid, otherwise may be undefined.
 */
DefInt SFBitsetAux3D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const DefInt region_length,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = 2 * region_length + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    DefInt index_min = region_length;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefInt iz = 0; iz <= region_length; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length - iz) * total_length  + region_length;
        for (DefInt iy = 0; iy <= region_length; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length + region_length;
            for (DefInt ix = 0; ix <= region_length; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        if (index_min > DefInt(ix)) {
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
                for (DefInt ix = 0; ix < region_length; ++ix) {
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
                            if (index_min > (ix + 1)) {
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
                    if (index_min > (iy)) {
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
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                index_min = 0;
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefInt iy = 0; iy < region_length; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                for (DefInt ix = 0; ix <= region_length; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > (ix)) {
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
                    for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                if (index_min > (ix + 1)) {
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
                        if (index_min > (iy + 1)) {
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
                if (index_min > (iz)) {
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
        for (DefInt iz = 0; iz < region_length; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length + iz + 1) * total_length + region_length;
            for (DefInt iy = 0; iy <= region_length; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length + region_length;
                for (DefInt ix = 0; ix <= region_length; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > (ix)) {
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
                    for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                if (index_min > (ix + 1)) {
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
                        if (index_min > (iy)) {
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
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    index_min = 0;
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefInt iy = 0; iy < region_length; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length + region_length;
                    for (DefInt ix = 0; ix <= region_length; ++ix) {
                        vec_index_x = vec_index_y - ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            if (index_min > (ix)) {
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
                        for (DefInt ix = 0; ix < region_length; ++ix) {
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
                                    if (index_min > (ix + 1)) {
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
                            if (index_min > (iy + 1)) {
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
                    if (index_min > (iz + 1)) {
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
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length_neg the length of searching region in negative directions.
 * @param[in] region_length_pos the length of searching region in positive directions.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 */
void SFBitsetAux3D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const DefInt region_length_neg, const DefInt region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length = region_length_neg + region_length_pos + 1;
    ptr_sfbitset_nodes->resize(total_length * total_length * total_length);
    ptr_sfbitset_nodes->assign(total_length * total_length * total_length, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefInt iz = 0; iz <= region_length_neg; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length_neg - iz) * total_length  + region_length_neg;
        for (DefInt iy = 0; iy <= region_length_neg; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length + region_length_neg;
            for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length + region_length_neg;
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefInt iy = 0; iy < region_length_pos; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length + region_length_neg;
                for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length + region_length_neg;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
            bool_not_z_max = false;
        }
    }
    if (bool_not_z_max) {
        for (DefInt iz = 0; iz < region_length_pos; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length_neg + iz + 1) * total_length + region_length_neg;
            for (DefInt iy = 0; iy <= region_length_neg; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length + region_length_neg;
                for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length + region_length_neg;
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefInt iy = 0; iy < region_length_pos; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length + region_length_neg;
                    for (DefInt ix = 0; ix <= region_length_neg; ++ix) {
                        vec_index_x = vec_index_y - ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length + region_length_neg;
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefInt ix = 0; ix < region_length_pos; ++ix) {
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
                    break;
                }
            }
        }
    }
}
/**
 * @brief function to find nodes nearby within a given distance.
 * @param[in] sfbitset_in the input node.
 * @param[in] region_length_neg lengths of searching region in negative directions.
 * @param[in] region_length_pos lengths of searching region in positive directions.
 * @param[in] periodic_min booleans indicating if the minimum domain boundaries are periodic.
 * @param[in] periodic_max booleans indicating if the maximum domain boundaries are periodic.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[out] ptr_sfbitset_nodes pointer to found nodes.
 */
void SFBitsetAux3D::FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
    const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level,
    const std::vector<DefSFBitset>& domain_max_n_level,
    std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const {
#ifdef DEBUG_CHECK_GRID
    if (domain_min_n_level.size() != 3 || domain_max_n_level.size() != 3) {
        LogManager::LogError("dimension of minimum (" + std::to_string(domain_min_n_level.size())
             + ") or maximum (" + std::to_string(domain_max_n_level.size())
             + ") boundary is not equal to space filling code dimension (3)");
    }
    if (region_length_neg.size() < 3 || region_length_pos.size() < 3) {
        LogManager::LogError("size of input searching length should not be less than the dimension");
    }
#endif  // DEBUG_CHECK_GRID
    DefAmrLUint total_length_x = region_length_neg[kXIndex] + region_length_pos[kXIndex] + 1;
    DefAmrLUint total_length_y = region_length_neg[kYIndex] + region_length_pos[kYIndex] + 1;
    DefAmrLUint total_length_z = region_length_neg[kZIndex] + region_length_pos[kZIndex] + 1;
    ptr_sfbitset_nodes->resize(total_length_x * total_length_y * total_length_z);
    ptr_sfbitset_nodes->assign(total_length_x * total_length_y * total_length_z, kInvalidSFbitset);
    DefSFBitset sfbitset_tmp_z = sfbitset_in, sfbitset_tmp_y, sfbitset_tmp_x;
    DefAmrLUint vec_index_x, vec_index_y, vec_index_z;
    bool bool_not_x_max, bool_not_y_max, bool_not_z_max;
    // negative z direction
    for (DefInt iz = 0; iz <= region_length_neg[kZIndex]; ++iz) {
        sfbitset_tmp_y = sfbitset_tmp_z;
        vec_index_z = (region_length_neg[kZIndex] - iz) * total_length_y  + region_length_neg[kYIndex];
        for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
            for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                vec_index_x = vec_index_y - ix;
                ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_min_n_level.at(kXIndex)) {
                    if (periodic_min.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_max_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                    } else {
                        break;
                    }
                }
                sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
            }
            sfbitset_tmp_x = sfbitset_tmp_y;
            vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
            bool_not_x_max = true;
            if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                if (periodic_max.at(kXIndex)) {
                    sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kXIndex);
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                } else {
                    bool_not_x_max = false;
                }
            }
            if (bool_not_x_max) {
                for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
                sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                    |domain_min_n_level.at(kYIndex);
                sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
            } else {
                bool_not_y_max = false;
            }
        }
        if (bool_not_y_max) {
            for (DefInt iy = 0; iy < region_length_pos[kYIndex]; ++iy) {
                sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length_x + region_length_neg[kXIndex];
                for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z + iy + 1) * total_length_x + region_length_neg[kXIndex];
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
            bool_not_z_max = false;
        }
    }
    if (bool_not_z_max) {
        for (DefInt iz = 0; iz < region_length_pos[kZIndex]; ++iz) {
            sfbitset_tmp_z = FindZPos(sfbitset_tmp_z);
            sfbitset_tmp_y = sfbitset_tmp_z;
            vec_index_z = (region_length_neg[kZIndex] + iz + 1) * total_length_y + region_length_neg[kYIndex];
            for (DefInt iy = 0; iy <= region_length_neg[kYIndex]; ++iy) {
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
                for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                    vec_index_x = vec_index_y - ix;
                    ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                    sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                }
                sfbitset_tmp_x = sfbitset_tmp_y;
                vec_index_y = (vec_index_z - iy) * total_length_x + region_length_neg[kXIndex];
                bool_not_x_max = true;
                if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                    == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                    if (periodic_max.at(kXIndex)) {
                        sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                            |domain_min_n_level.at(kXIndex);
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    } else {
                        bool_not_x_max = false;
                    }
                }
                if (bool_not_x_max) {
                    for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
                    sfbitset_tmp_y = (sfbitset_tmp_y&k0SFBitsetTakeYRef_.at(kRefOthers_))
                        |domain_min_n_level.at(kYIndex);
                    sfbitset_tmp_y = FindYNeg(sfbitset_tmp_y);
                } else {
                    bool_not_y_max = false;
                }
            }
            if (bool_not_y_max) {
                for (DefInt iy = 0; iy < region_length_pos[kYIndex]; ++iy) {
                    sfbitset_tmp_y = FindYPos(sfbitset_tmp_y);
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length_x + region_length_neg[kXIndex];
                    for (DefInt ix = 0; ix <= region_length_neg[kXIndex]; ++ix) {
                        vec_index_x = vec_index_y - ix;
                        ptr_sfbitset_nodes->at(vec_index_x) = sfbitset_tmp_x;
                        if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_min_n_level.at(kXIndex)) {
                        if (periodic_min.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_max_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXPos(sfbitset_tmp_x);
                        } else {
                            break;
                        }
                    }
                        sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                    }
                    sfbitset_tmp_x = sfbitset_tmp_y;
                    vec_index_y = (vec_index_z + iy + 1) * total_length_x + region_length_neg[kXIndex];
                    bool_not_x_max = true;
                    if ((sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefCurrent_))
                        == domain_max_n_level.at(kXIndex)) {  // if the start node at max x boundary
                        if (periodic_max.at(kXIndex)) {
                            sfbitset_tmp_x = (sfbitset_tmp_x&k0SFBitsetTakeXRef_.at(kRefOthers_))
                                |domain_min_n_level.at(kXIndex);
                            sfbitset_tmp_x = FindXNeg(sfbitset_tmp_x);
                        } else {
                            bool_not_x_max = false;
                        }
                    }
                    if (bool_not_x_max) {
                        for (DefInt ix = 0; ix < region_length_pos[kXIndex]; ++ix) {
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
                    break;
                }
            }
        }
    }
}
/**
 * @brief function to check if node at current level should exist based on lower level nodes.
 * @param[in] sfbitset_in input node at current level.
 * @param[in] exist_nodes_lower_level existing nodes at one lower level.
 * @return true if the node at current level should exist.
 */
bool SFBitsetAux3D::CheckExistenceCurrentLevel(
    const DefSFBitset& sfbitset_in, const DefMap<DefInt>& exist_nodes_lower_level) const {
    DefSFBitset sfbitset_in_lower = SFBitsetToOneLowerLevel(sfbitset_in);
    if (exist_nodes_lower_level.find(sfbitset_in_lower) == exist_nodes_lower_level.end()) {
        return false;
    }
    DefSFBitset current_bit = sfbitset_in&k0SfBitsetCurrentLevelBits_;
    bool bool_x = (current_bit&k0SFBitsetTakeXRef_.at(kRefCurrent_)) != 0,
        bool_y = (current_bit&k0SFBitsetTakeYRef_.at(kRefCurrent_)) != 0,
        bool_z = (current_bit&k0SFBitsetTakeZRef_.at(kRefCurrent_)) != 0;
    if (current_bit != 0) {
        if (bool_x) {
            DefSFBitset sfbitset0 = FindXPos(sfbitset_in_lower);
            if (exist_nodes_lower_level.find(sfbitset0)
                == exist_nodes_lower_level.end()) {
                return false;
            }
            DefSFBitset sfbitset1 = FindYPos(sfbitset0);
            if (bool_y) {
                if (exist_nodes_lower_level.find(sfbitset1)
                    == exist_nodes_lower_level.end()) {
                    return false;
                }
                if (bool_z) {
                    if (exist_nodes_lower_level.find(FindZPos(sfbitset1))
                        == exist_nodes_lower_level.end()) {
                        return false;
                    }
                }
            }
            if (bool_z) {
                if (exist_nodes_lower_level.find(FindZPos(sfbitset0))
                    == exist_nodes_lower_level.end()) {
                    return false;
                }
            }
        }
        if (bool_y) {
            DefSFBitset sfbitset0 = FindYPos(sfbitset_in_lower);
            if (exist_nodes_lower_level.find(sfbitset0)
                == exist_nodes_lower_level.end()) {
                return false;
            }
            if (bool_z) {
                if (exist_nodes_lower_level.find(FindZPos(sfbitset0))
                    == exist_nodes_lower_level.end()) {
                    return false;
                }
            }
        }
        if (bool_z) {
            if (exist_nodes_lower_level.find(FindYPos(sfbitset_in_lower))
                == exist_nodes_lower_level.end()) {
                return false;
            }
        }
    }
    return true;
}
/**
 * @brief function to check if the node is on a given boundary.
 * @param[in] i_dir the given direction.
 * @param[in] sfbitset_in space filling code of a node.
 * @param[in] sfbitset_boundary space filling code of the boundary in the given direction.
 * @return true if the node is on the given boundary.
 */
bool SFBitsetAux3D::CheckIfOnGivenBoundary(const DefInt i_dir,
    const DefSFBitset& sfbitset_in, const DefSFBitset& sfbitset_boundary) const {
    switch (i_dir) {
    case kXIndex:
        if ((sfbitset_in&k0SFBitsetTakeXRef_.at(kRefCurrent_)) == sfbitset_boundary) {
            return true;
        } else {
            return false;
        }
    case kYIndex:
        if ((sfbitset_in&k0SFBitsetTakeYRef_.at(kRefCurrent_)) == sfbitset_boundary) {
            return true;
        } else {
            return false;
        }
    case kZIndex:
        if ((sfbitset_in&k0SFBitsetTakeZRef_.at(kRefCurrent_)) == sfbitset_boundary) {
            return true;
        } else {
            return false;
        }
    default:
        return false;
    }
}
/**
 * @brief function to calculate spacing fill code of minimum indices minus 1 at a given level.
 * @param[in] i_level the given refinement level
 * @param[in] indices_min minimum indicies of the computational domain 
 * @param[out] ptr_min_m1_bitsets a pointer to minimum indices minus 1
 * @throws ErrorType if the size of min_m1_bitsets is not 3
 */
void SFBitsetAux3D::GetMinM1AtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_min,
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
void SFBitsetAux3D::GetMaxP1AtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_max,
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
void SFBitsetAux3D::GetMinAtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_min,
    std::vector<DefSFBitset>* const ptr_min_bitsets) const {
    ptr_min_bitsets->resize(3);
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
void SFBitsetAux3D::GetMaxAtGivenLevel(const DefInt i_level,
    std::vector<DefAmrLUint> indices_max,
    std::vector<DefSFBitset>* const ptr_max_bitsets) const {
    ptr_max_bitsets->resize(3);
    ptr_max_bitsets->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({indices_max[kXIndex], 0, 0}));
    ptr_max_bitsets->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, indices_max[kYIndex], 0}));
    ptr_max_bitsets->at(kZIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetEncoding({0, 0, indices_max[kZIndex]}));
}
/**
* @brief   function to find space filling code of nodes on a cell (3D) which may across two refinement levels
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  map_node_exist nodes at fine grid with i_level - 1 spacing filling code
* @param[in]  map_node_exist_coarse  nodes at coarse grid with i_level - 1 spacing filling code
* @param[out] ptr_sfbitsets nodes on the cell 
* @return  if true indicates the given node belongs to a cell
* @note ptr_sfbitsets[0]:(0, 0, 0); ptr_sfbitsets[1]:(+x, 0, 0);
*       ptr_sfbitsets[2]:(0, +y, 0); ptr_sfbitsets[3]:(+x, +y, 0);
*       ptr_sfbitsets[4]:(0, 0, +z);  ptr_sfbitsets[5]:(+x, 0, +z);
*       ptr_sfbitsets[6]:(0, +y, +z); ptr_sfbitsets[7]: (+x, +y, +z).
*/
bool SFBitsetAux3D::SFBitsetBelongToOneCellAcrossTwoLevels(const DefSFBitset& sfbitset_in,
    const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_node_exist_coarse,
    std::array<std::pair<DefSFBitset, DefInt>, 8>* const ptr_sfbitsets) const {
    if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(sfbitset_in) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(sfbitset_in) == 0) {
            return false;
        }
        ptr_sfbitsets->at(0).second = 2;
    } else {
        ptr_sfbitsets->at(0).second = 1;
    }
    ptr_sfbitsets->at(0).first = sfbitset_in;
    // (+x, 0, 0)
    ptr_sfbitsets->at(1).first = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(1).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(1).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(1).second = 2;
    } else {
        ptr_sfbitsets->at(1).second = 1;
    }
    // (+x, +y, 0)
    ptr_sfbitsets->at(3).first = FindYPos(ptr_sfbitsets->at(1).first);
    if (map_node_exist.find(ptr_sfbitsets->at(3).first)  == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(3).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(3).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(3).second = 2;
    } else {
        ptr_sfbitsets->at(3).second = 1;
    }
    // (0, +y, 0)
    ptr_sfbitsets->at(2).first = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2).first)  == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(2).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(2).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(2).second = 2;
    } else {
        ptr_sfbitsets->at(2).second = 1;
    }
    // (0, 0, +z)
    ptr_sfbitsets->at(4).first = FindZPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(4).first)  == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(4).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(4).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(4).second = 2;
    } else {
        ptr_sfbitsets->at(4).second = 1;
    }
    // (+x, 0, +z)
    ptr_sfbitsets->at(5).first = FindXPos(ptr_sfbitsets->at(4).first);
    if (map_node_exist.find(ptr_sfbitsets->at(5).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(5).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(5).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(5).second = 2;
    } else {
        ptr_sfbitsets->at(5).second = 1;
    }
    // (+x, +y, +z)
    ptr_sfbitsets->at(7).first = FindYPos(ptr_sfbitsets->at(5).first);
    if (map_node_exist.find(ptr_sfbitsets->at(7).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(7).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(7).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(7).second = 2;
    } else {
        ptr_sfbitsets->at(7).second = 1;
    }
    // (0, +y, +z)
    ptr_sfbitsets->at(6).first = FindYPos(ptr_sfbitsets->at(4).first);
    if (map_node_exist.find(ptr_sfbitsets->at(6).first) == map_node_exist.end()) {
        if (map_node_exist_coarse.find(ptr_sfbitsets->at(6).first) == map_node_exist_coarse.end()
            || map_node_exist_coarse.at(ptr_sfbitsets->at(6).first) == 0) {
            return false;
        }
        ptr_sfbitsets->at(6).second = 2;
    } else {
        ptr_sfbitsets->at(6).second = 1;
    }
    return true;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
