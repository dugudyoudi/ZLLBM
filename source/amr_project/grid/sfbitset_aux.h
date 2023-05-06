//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file sfbitset_aux.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define class and delare functions for space fill code.
*/

#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
#include <array>
#include <vector>
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
class  SFBitsetAuxInterface {
 public:
    // .at(kRefCurrent) takes the given direction while
    // .at(kRefOthers) takes directions other than the given one
    static constexpr DefUint kRefOthers_ = 0, kRefCurrent_ = 1;
    DefSFBitset k0SfBitsetCurrentLevelBits_ = 0;
    /**< reference bitset used to take digitals at the center of
    a cell of one level lower.*/
    std::vector<DefReal> k0SpaceBackground_;
    virtual void SFBitsetInitial() = 0;
    virtual DefSFBitset SFBitsetEncodingCoordi(
        const std::vector<DefReal>& grid_space,
        const std::vector<DefReal>& coordi) const = 0;
    // compared to inline function, this virtual is costly, avoiding using this if possible
    virtual DefSFBitset SFBitsetToNLowerLevelVir(
        const DefSizet n_level, const DefSFBitset& morton_in) const = 0;
    virtual ~SFBitsetAuxInterface() {}
};
#ifndef  DEBUG_DISABLE_2D_FUNCTION
class  SFBitsetAux2D : public SFBitsetAuxInterface {
 public:
    static constexpr DefUint  kNodeIndexX0Y0_ = 0, kNodeIndexXnY0_ = 1,
        kNodeIndexXpY0_ = 2, kNodeIndexX0Yn_ = 3, kNodeIndexX0Yp_ = 4,
        kNodeIndexXnYn_ = 5, kNodeIndexXnYp_ = 6, kNodeIndexXpYn_ = 7,
        kNodeIndexXpYp_ = 8;
    ///< indices of node and its 8 neighbours
    static std::array<DefSFBitset, 2>
        k0SFBitsetTakeXRef_, k0SFBitsetTakeYRef_;
    /**< reference bitset used to take digitals at a given direction
          by bool operator.*/
    std::array<DefSFBitset, 2> SFBitsetMin_, SFBitsetMax_;
    std::array<DefSFBitset, 2> k0SFBitsetGlobleMin_, k0SFBitsetGlobleMax_;
    /**< bitset corresponding to the minimum and maximum
     coordinate in each direction*/

    /* space filling code related functions: when use code other than morton
    code, need to rewrite the following functions and those in
    (morton_aux.cpp). Virtual functions are not used here in oder to
    avoid extra expends. Functions will be called frequently are declared
    as inline functions.*/
    inline DefSFBitset SFBitsetToOneLowerLevel(
        const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToOneHigherLevel(
        const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNLowerLevel(
        const DefSizet n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNHigherLevel(
        const DefSizet n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset FindXNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindXPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYPos(const DefSFBitset& bitset_in) const;

    void SFBitsetComputeIndices(const DefSFBitset& bitset_in,
        std::array<DefLUint, 2>* const ptr_indiex) const;
    void SFBitsetComputeCoordinate(
        const DefSFBitset& bitset_in,
        const std::array<DefReal, 2>& grid_space,
        std::array<DefReal, 2>* const ptr_coordi) const;
    DefSFBitset SFBitsetEncoding(
        const std::array<DefLUint, 2>& coordi_index) const;

    DefSFBitset SFBitsetBitsForRefinement(const DefSizet i_level) const;
    void SFBitsetMinAndMaxCoordinates(const DefSizet max_level,
        const std::array<DefLUint, 2>& indices_min,
        const std::array<DefLUint, 2>& indices_max);
    void SFBitsetMinAndMaxGlobal(
        const std::array<DefLUint, 2>& indices_min,
        const std::array<DefLUint, 2>& indices_max);
    void SFBitsetNotOnDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::array<DefSFBitset, 2>& sfbitset_min,
        const std::array<DefSFBitset, 2>& sfbitset_max,
        std::array<bool, 2>* const ptr_bool_not_at_boundary_neg,
        std::array<bool, 2>* const ptr_bool_not_at_boundary_pos) const;
    void SFBitsetFindAllNeighbours(
            const DefSFBitset& sfbitset_center,
            std::array<DefSFBitset, 9>* const ptr_bitset_neighbour) const;
    void SFBitsetFindCellNeighbours(
        const DefSFBitset& sfbitset_corner,
        std::array<DefSFBitset, 4>* const ptr_bitset_neighbour) const;
    void SFBitsetHigherLevelInACell(
        const DefSizet level_diff,
        const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_vec_bitset_neighbour) const;
    template<typename DataType>
    bool SFBitsetBelongToOneCell(const DefSFBitset& sfbitset_in,
        const DefMap<DataType>& map_node_exist,
        std::array<DefSFBitset, 4>* const ptr_sfbitsets) const;

    // virtual functions
    void SFBitsetInitial() override;
    DefSFBitset SFBitsetEncodingCoordi(
        const std::vector<DefReal>& grid_space,
        const std::vector<DefReal>& coordi) const override;
    DefSFBitset  SFBitsetToNLowerLevelVir(const DefSizet n_level, const DefSFBitset& morton_in) const final {
        return SFBitsetToNLowerLevel(n_level, morton_in);
    }

    SFBitsetAux2D() { SFBitsetInitial(); }
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
class  SFBitsetAux3D : public SFBitsetAuxInterface {
 public:
    // p is positive and n is negative
    static constexpr DefUint  kNodeIndexX0Y0Z0_ = 0, kNodeIndexXnY0Z0_ = 1,
        kNodeIndexXpY0Z0_ = 2, kNodeIndexX0YnZ0_ = 3, kNodeIndexX0YpZ0_ = 4,
        kNodeIndexX0Y0Zn_ = 9, kNodeIndexX0Y0Zp_ = 10, kNodeIndexXnYnZ0_ = 5,
        kNodeIndexXnYpZ0_ = 6, kNodeIndexXpYnZ0_ = 7, kNodeIndexXpYpZ0_ = 8,
        kNodeIndexXnY0Zn_ = 11, kNodeIndexXnY0Zp_ = 12, kNodeIndexXpY0Zn_ = 13,
        kNodeIndexXpY0Zp_ = 14, kNodeIndexX0YnZn_ = 15, kNodeIndexX0YnZp_ = 16,
        kNodeIndexX0YpZn_ = 17, kNodeIndexX0YpZp_ = 18, kNodeIndexXpYpZp_ = 19,
        kNodeIndexXpYnZp_ = 20, kNodeIndexXpYpZn_ = 21, kNodeIndexXpYnZn_ = 22,
        kNodeIndexXnYpZp_ = 23, kNodeIndexXnYnZp_ = 24, kNodeIndexXnYpZn_ = 25,
        kNodeIndexXnYnZn_ = 26;  ///< indices of node and  its 26 neighbours
    static std::array<DefSFBitset, 2>
        k0SFBitsetTakeXRef_, k0SFBitsetTakeYRef_, k0SFBitsetTakeZRef_;
    /**< reference bitset used to take digitals at a given direction
          by bool operator.*/
    std::array<DefSFBitset, 3> SFBitsetMin_, SFBitsetMax_;
    std::array<DefSFBitset, 3> k0SFBitsetGlobleMin_, k0SFBitsetGlobleMax_;
    /**< bitset corresponding to the minimum and maximum
     coordinate in each direction*/
    inline DefSFBitset SFBitsetToOneLowerLevel(
        const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToOneHigherLevel(
        const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNLowerLevel(
        const DefSizet n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNHigherLevel(
        const DefSizet n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset FindXNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindXPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindZNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindZPos(const DefSFBitset& bitset_in) const;

    DefSFBitset SFBitsetEncoding(
        const std::array<DefLUint, 3>& coordi_index) const;
    void SFBitsetComputeIndices(const DefSFBitset& bitset_in,
        std::array<DefLUint, 3>* const ptr_indices) const;
    void SFBitsetComputeCoordinate(const DefSFBitset& bitset_in,
        const std::array<DefReal, 3>& grid_space,
        std::array<DefReal, 3>* const ptr_coordi) const;

    DefSFBitset SFBitsetBitsForRefinement(const DefSizet i_level) const;
    void SFBitsetMinAndMaxCoordinates(const DefSizet max_level,
        const std::array<DefLUint, 3>& indices_min,
        const std::array<DefLUint, 3>& indices_max);
    void SFBitsetMinAndMaxGlobal(
        const std::array<DefLUint, 3>& indices_min,
        const std::array<DefLUint, 3>& indices_max);
    void SFBitsetNotOnDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::array<DefSFBitset, 3>& sfbitset_min,
        const std::array<DefSFBitset, 3>& sfbitset_max,
        std::array<bool, 3>* const ptr_bool_not_at_boundary_neg,
        std::array<bool, 3>* const ptr_bool_not_at_boundary_pos) const;
    void SFBitsetFindAllNeighbours(
        const DefSFBitset& sfbitset_center,
        std::array<DefSFBitset, 27>* const ptr_bitset_neighbour) const;
    void SFBitsetFindCellNeighbours(
        const DefSFBitset& sfbitset_corner,
        std::array<DefSFBitset, 8>* const ptr_vec_bitset_neighbour) const;
    void SFBitsetHigherLevelInACell(
        const DefSizet level_diff,
        const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_ptr_sfbitsets_higher_level) const;
    template<typename DataType>
    bool SFBitsetBelongToOneCell(const DefSFBitset& sfbitset_in,
        const DefMap<DataType>& map_node_exist,
        std::array<DefSFBitset, 8>* const ptr_sfbitsets) const;

    // virtual functions
    void SFBitsetInitial() override;
    DefSFBitset SFBitsetEncodingCoordi(
        const std::vector<DefReal>& grid_space,
        const std::vector<DefReal>& coordi) const override;
        DefSFBitset  SFBitsetToNLowerLevelVir(const DefSizet n_level, const DefSFBitset& morton_in) const final {
        return SFBitsetToNLowerLevel(n_level, morton_in);
    }

    SFBitsetAux3D() { SFBitsetInitial(); }
};
#endif  // DEBUG_DISABLE_3D_FUNCTION
}  // end namespace amrproject
}  // end namespace rootproject

// include inline incomplement of morton code function
#include "grid/morton_aux.hpp"

// functions using space filling code manipulations
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to find space filling code of nodes on a cell (2D)
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  map_node_exist grid instance containing nodes at the same level
* @param[out] ptr_sfbitsets of nodes on the cell 
* @return  a bool value if true indicates
           the given node is at the origin of a cell
* @note ptr_sfbitsets[0]:(0, 0); ptr_sfbitsets[1]:(+x, 0);
*       ptr_sfbitsets[2]:(0, +y); ptr_sfbitsets[3]:(+x, +y);
*/
template<typename DataType>
bool SFBitsetAux2D::SFBitsetBelongToOneCell(
    const DefSFBitset& sfbitset_in, const DefMap<DataType>& map_node_exist,
    std::array<DefSFBitset, 4>* const ptr_sfbitsets)  const {

    bool bool_nodes_can_construct_a_cell = true;

    if (map_node_exist.find(sfbitset_in)
        == map_node_exist.end()) {
        return false;
    }
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1))
        == map_node_exist.end()) {
        return false;
    }
    // (+x, +y)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    if (map_node_exist.find(ptr_sfbitsets->at(3))
        == map_node_exist.end()) {
        return false;
    }
    // (0, +y)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2))
        == map_node_exist.end()) {
        return false;
    }
    return true;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to find space filling code constucting a cell (3D)
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  map_node_exist grid instance containing nodes at the same level
* @param[out] ptr_sfbitsets of nodes on the cell
* @return  a bool value if true indicates
           the given node is at the origin of a cell
* @note ptr_sfbitsets[0]:(0, 0, 0); ptr_sfbitsets[1]:(+x, 0, 0);
*       ptr_sfbitsets[2]:(0, +y, 0); ptr_sfbitsets[3]:(+x, +y, 0);
*       ptr_sfbitsets[4]:(0, 0, +z);  ptr_sfbitsets[5]:(+x, 0, +z);
*       ptr_sfbitsets[6]:(0, +y, +z); ptr_sfbitsets[7]: (+x, +y, +z).
*/
template<typename DataType>
bool SFBitsetAux3D::SFBitsetBelongToOneCell(
    const DefSFBitset& sfbitset_in, const DefMap<DataType>& map_node_exist,
    std::array<DefSFBitset, 8>* const ptr_sfbitsets) const {

    if (map_node_exist.find(sfbitset_in)
        == map_node_exist.end()) {
        return false;
    }
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1))
        == map_node_exist.end()) {
        return false;
    }
    // (+x, +y, 0)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    if (map_node_exist.find(ptr_sfbitsets->at(3))
        == map_node_exist.end()) {
        return false;
    }
    // (0, +y, 0)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2))
        == map_node_exist.end()) {
        return false;
    }
    // (0, 0, +z)
    ptr_sfbitsets->at(4) = FindZPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(4))
        == map_node_exist.end()) {
        return false;
    }
    // (+x, 0, +z)
    ptr_sfbitsets->at(5) = FindXPos(ptr_sfbitsets->at(4));
    if (map_node_exist.find(ptr_sfbitsets->at(5))
        == map_node_exist.end()) {
        return false;
    }
    // (+x, +y, +z)
    ptr_sfbitsets->at(7) = FindYPos(ptr_sfbitsets->at(5));
    if (map_node_exist.find(ptr_sfbitsets->at(7))
        == map_node_exist.end()) {
        return false;
    }
    // (0, +y, +z)
    ptr_sfbitsets->at(6) = FindYPos(ptr_sfbitsets->at(4));
    if (map_node_exist.find(ptr_sfbitsets->at(6))
        == map_node_exist.end()) {
        return false;
    }
    return true;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
