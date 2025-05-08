//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file sfbitset_aux.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define class and functions for space fill code.
*/

#ifndef SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
#define SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
#include <array>
#include <cstring>
#include <vector>
#include <utility>
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
/**
* @class SFBitsetAuxInterface
* @brief class used to handle functions for space filling curves
*/
class SFBitsetAuxInterface {
 public:
    static constexpr DefSFBitset kInvalidSFbitset = ~0;
    /**< reference bitset used to take digitals at a given direction
          by bool operator.*/
    // for k0SFBitsetTakeXRef_, k0SFBitsetTakeYRef_ and k0SFBitsetTakeZRef_
    // .at(kRefCurrent) takes the given direction (x, y, or z) while
    // .at(kRefOthers) takes directions other than the given one (y z, x z, or x y)
    static constexpr DefInt kRefOthers_ = 0, kRefCurrent_ = 1;
    DefSFBitset k0SfBitsetCurrentLevelBits_ = 0;
    /**< reference bitset used to take digitals at the center of
    a cell of one level lower.*/

    std::vector<DefReal> GetBackgroundGridSpacing() const {
        return k0SpaceBackground_;
    }

    inline DefSFCodeToUint SFBitsetToSFCode(const DefSFBitset& sfbitset) const { return sfbitset.to_ullong(); }

    virtual void SFBitsetInitial() = 0;
    virtual DefSFBitset SFBitsetEncodingCoordi(
        const std::vector<DefReal>& grid_space, const std::vector<DefReal>& coordi) const = 0;
    virtual void SFBitsetFindEdgeNode(const DefSFBitset& morton_in,
        std::vector<DefSFBitset>* const ptr_vec_edge_nodes) const = 0;
    virtual DefSFBitset FindNodeInNeg(const DefInt dir, const DefSFBitset& sfbitset_in) const = 0;
    virtual DefSFBitset FindNodeInPos(const DefInt dir, const DefSFBitset& sfbitset_in) const = 0;
    virtual void FindNodesInReginOfGivenLength(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const = 0;
    virtual DefInt FindNodesInPeriodicRegionCorner(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const = 0;
    virtual DefInt FindNodesInPeriodicRegionCornerOverlap(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const = 0;
    virtual DefInt FindNodesInPeriodicRegionCenterOverlap(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const = 0;
    virtual DefInt FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const = 0;
    virtual void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const = 0;
    virtual void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length_neg, const DefInt region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const = 0;
    virtual DefSFBitset SetIindexinGivenDirection(const DefInt i_dir, const DefSFBitset& sfbitset_in,
        const DefSFBitset& sfbitset_target) const = 0;
    virtual void SFBitsetHigherLevelInACell(
        const DefInt level_diff, const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_vec_bitset_neighbour) const = 0;

    virtual void FindNeighboringCoarseFromFine(const DefSFBitset& sfbitset_fine,
        std::vector<DefSFBitset>* const ptr_vec_coarse) const = 0;
    virtual bool CheckExistenceCurrentLevel(
        const DefSFBitset& sfbitset_in, const DefMap<DefInt>& exist_nodes) const = 0;
    virtual bool CheckIfOnGivenBoundary(const DefInt i_dir,
        const DefSFBitset& sfbitset_in, const DefSFBitset& sfbitset_boundary) const = 0;

    // compared to inline function, this virtual is costly, avoiding using this if possible
    virtual DefSFBitset SFBitsetToNLowerLevelVir(
        const DefInt n_level, const DefSFBitset& sfbitset_in) const = 0;
    virtual DefSFBitset SFBitsetToNHigherLevelVir(
        const DefInt n_level, const DefSFBitset& sfbitset_in) const = 0;
    virtual DefInt SFBitsetCoincideLevelVir(const DefSFBitset& sfbitset_in) const = 0;
    virtual void SFBitsetFindAllNeighborsVir(const DefSFBitset& sfbitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const = 0;
    virtual void SFBitsetFindAllBondedNeighborsVir(const DefSFBitset& sfbitset_in,
        const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const = 0;
    virtual void SFBitsetComputeCoordinateVir(const DefSFBitset& bitset_in,
        const std::vector<DefReal>& grid_space, std::vector<DefReal>* const ptr_coordi) const = 0;
    virtual void SFBitsetFindCellNeighborsVir(const DefSFBitset& sfbitset_corner,
            std::vector<DefSFBitset>* const ptr_bitset_neighbour) const = 0;

    virtual std::vector<DefAmrLUint> GetMinBackgroundIndices() const = 0;
    virtual std::vector<DefAmrLUint> GetMaxBackgroundIndices() const = 0;
    virtual void GetMinM1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const = 0;
    virtual void GetMaxP1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_p1_bitset) const = 0;
    virtual void GetMinAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_bitsets) const = 0;
    virtual void GetMaxAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_bitset) const = 0;

    virtual ~SFBitsetAuxInterface() {}

    // set and get protected members
    virtual void SetSpaceBackground(const std::vector<DefReal>& space_background) = 0;
    std::array<DefSFBitset, 2> GetTakeXRef() const { return k0SFBitsetTakeXRef_; }
    std::array<DefSFBitset, 2> GetTakeYRef() const { return k0SFBitsetTakeYRef_; }
    std::array<DefSFBitset, 2> GetTakeZRef() const { return k0SFBitsetTakeZRef_; }

    // mpi related
    virtual void GetNLevelCorrespondingOnes(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_last_ones) const = 0;

 protected:
    std::vector<DefReal> k0SpaceBackground_;
    std::array<DefSFBitset, 2> k0SFBitsetTakeXRef_, k0SFBitsetTakeYRef_, k0SFBitsetTakeZRef_;
    const DefInt max_reset_code_ = 1000;
    ///<  maximum iteration for reseting indices
};
#ifndef  DEBUG_DISABLE_2D_FUNCTION
class  SFBitsetAux2D : public SFBitsetAuxInterface {
 public:
    static constexpr DefInt  kNodeIndexX0Y0_ = 0, kNodeIndexXnY0_ = 1,
     kNodeIndexXpY0_ = 2, kNodeIndexX0Yn_ = 3, kNodeIndexX0Yp_ = 4,
     kNodeIndexXnYn_ = 5, kNodeIndexXnYp_ = 6, kNodeIndexXpYn_ = 7,
     kNodeIndexXpYp_ = 8;
    ///< indices of node and its 8 neighbors
    std::array<DefSFBitset, 2> SFBitsetMin_ = {~DefSFCodeToUint(0), ~DefSFCodeToUint(0)}, SFBitsetMax_ = {0, 0};
    /**< bitset corresponding to the minimum and maximum
     indices of the computational domain in each direction*/

    /* give an offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (k0MinIndexOfBackgroundNode_ * kDomainDx),
    and the default value is 1. */
    std::array<DefAmrLUint, 2> k0MinIndexOfBackgroundNode_ = {1, 1};
    ///< the minimum index of background nodes in each direction*/
    std::array<DefAmrLUint, 2> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0MinIndexOfBackgroundNode_ is included in k0MaxIndexOfBackgroundNode_
    // k0DomainSize_ is k0DomainDx_ * (k0MaxIndexOfBackgroundNode_ - k0MinIndexOfBackgroundNode_)
    std::array<DefReal, 2> k0DomainSize_{};  ///< domain size
    std::array<DefReal, 2> k0DomainDx_{};  ///< grid space
    std::array<DefReal, 2> k0RealMin_{};  ///< k0MinIndexOfBackgroundNode_ * kDomainDx

    /* space filling code related functions: when use code other than morton
    code, need to rewrite the following functions and those in
    (morton_aux.cpp). Virtual functions are not used here in oder to
    avoid extra expends. Functions will be called frequently are declared
    as inline functions.*/
    inline DefSFBitset SFBitsetToOneLowerLevel(const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToOneHigherLevel(const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNLowerLevel(const DefInt n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNHigherLevel(const DefInt n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset FindXNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindXPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYPos(const DefSFBitset& bitset_in) const;
    inline DefInt SFBitsetCoincideLevel(const DefSFBitset& bitset_in) const;

    void SFBitsetComputeIndices(const DefSFBitset& bitset_in,
        std::array<DefAmrLUint, 2>* const ptr_indiex) const;
    void SFBitsetComputeCoordinate(const DefSFBitset& bitset_in,
        const std::array<DefReal, 2>& grid_space, std::array<DefReal, 2>* const ptr_coordi) const;
    DefSFBitset SFBitsetEncoding(const std::array<DefAmrLUint, 2>& coordi_index) const;

    DefSFBitset SFBitsetBitsForRefinement(const DefInt i_level) const;
    void SFBitsetSetMinAndMaxBounds(
        const std::array<DefAmrLUint, 2>& indices_min, const std::array<DefAmrLUint, 2>& indices_max);
    void SFBitsetNotOnDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::array<DefSFBitset, 2>& sfbitset_min, const std::array<DefSFBitset, 2>& sfbitset_max,
        std::array<bool, 2>* const ptr_bool_not_at_boundary_neg,
        std::array<bool, 2>* const ptr_bool_not_at_boundary_pos) const;
    void SFBitsetFindAllNeighbors(const DefSFBitset& sfbitset_center,
        std::array<DefSFBitset, 9>* const ptr_bitset_neighbour) const;
    void SFBitsetFindCellNeighbors(const DefSFBitset& sfbitset_corner,
        std::array<DefSFBitset, 4>* const ptr_bitset_neighbour) const;
    template<typename DataType>
    bool SFBitsetBelongToOneCell(const DefSFBitset& sfbitset_in,
        const DefMap<DataType>& map_node_exist,
        std::array<DefSFBitset, 4>* const ptr_sfbitsets) const;
    bool SFBitsetBelongToOneCellAcrossTwoLevels(const DefSFBitset& sfbitset_in,
        const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_node_exist_coarse,
        std::array<std::pair<DefSFBitset, DefInt>, 4>* const ptr_sfbitsets) const;
    int ResetIndicesExceedingDomain(const std::array<DefAmrLUint, 2>& domain_min_indices,
        const std::array<DefAmrLUint, 2>& domain_max_indices,
        DefSFCodeToUint* const ptr_i_code, DefSFBitset* ptr_bitset_tmp) const;

    // virtual functions
    void SFBitsetInitial() override;
    DefSFBitset SFBitsetEncodingCoordi(const std::vector<DefReal>& grid_space,
        const std::vector<DefReal>& coordi) const override;
     void SFBitsetFindEdgeNode(const DefSFBitset& morton_in,
        std::vector<DefSFBitset>* const ptr_vec_edge_nodes) const final {
        ptr_vec_edge_nodes->resize(4);
        ptr_vec_edge_nodes->at(0) = FindXNeg(morton_in);
        ptr_vec_edge_nodes->at(1) = FindXPos(morton_in);
        ptr_vec_edge_nodes->at(2) = FindYNeg(morton_in);
        ptr_vec_edge_nodes->at(3) = FindYPos(morton_in);
    }
    DefSFBitset FindNodeInNeg(const DefInt dir, const DefSFBitset& sfbitset_in) const override {
        if (dir == kXIndex) {
            return FindXNeg(sfbitset_in);
        } else if (dir == kYIndex) {
            return FindYNeg(sfbitset_in);
        } else {
            return sfbitset_in;
        }
    }
    DefSFBitset FindNodeInPos(const DefInt dir, const DefSFBitset& sfbitset_in) const override {
        if (dir == kXIndex) {
            return FindXPos(sfbitset_in);
        } else if (dir == kYIndex) {
            return FindYPos(sfbitset_in);
        } else {
            return sfbitset_in;
        }
    }
    void SFBitsetComputeCoordinateVir(const DefSFBitset& bitset_in,
        const std::vector<DefReal>& grid_space, std::vector<DefReal>* const ptr_coordi) const final {
        std::array<DefReal, 2> arr_grid_space = {grid_space[kXIndex], grid_space[kYIndex]}, arr_coordi;
        SFBitsetComputeCoordinate(bitset_in, arr_grid_space, &arr_coordi);
        ptr_coordi->resize(2);
        ptr_coordi->at(kXIndex) = arr_coordi[kXIndex];
        ptr_coordi->at(kYIndex) = arr_coordi[kYIndex];
    }
    void SFBitsetFindCellNeighborsVir(const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_bitset_neighbour) const final {
         std::array<DefSFBitset, 4> array_neighbors;
        SFBitsetFindCellNeighbors(sfbitset_corner, &array_neighbors);
        ptr_bitset_neighbour->resize(4);
        memcpy(ptr_bitset_neighbour->data(), array_neighbors.data(),
            4 * sizeof(DefSFBitset));
    }
    void FindNodesInReginOfGivenLength(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const final;
    DefInt FindNodesInPeriodicRegionCorner(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    DefInt FindNodesInPeriodicRegionCornerOverlap(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const final;
    DefInt FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length_neg, const DefInt region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    DefInt FindNodesInPeriodicRegionCenterOverlap(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const final;
    DefSFBitset SetIindexinGivenDirection(const DefInt i_dir, const DefSFBitset& sfbitset_in,
        const DefSFBitset& sfbitset_target) const final;
    void SFBitsetHigherLevelInACell(
            const DefInt level_diff, const DefSFBitset& sfbitset_corner,
            std::vector<DefSFBitset>* const ptr_vec_bitset_neighbour) const final;
    void SFBitsetFindAllBondedNeighborsVir(const DefSFBitset& bitset_in,
        const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const final;

    bool CheckExistenceCurrentLevel(
        const DefSFBitset& sfbitset_in, const DefMap<DefInt>& exist_nodes) const final;
    bool CheckIfOnGivenBoundary(const DefInt i_dir,
        const DefSFBitset& sfbitset_in, const DefSFBitset& sfbitset_boundary) const final;

    DefSFBitset SFBitsetToNLowerLevelVir(const DefInt n_level,
        const DefSFBitset& morton_in) const final {
        return SFBitsetToNLowerLevel(n_level, morton_in);
    }
    DefSFBitset  SFBitsetToNHigherLevelVir(const DefInt n_level, const DefSFBitset& morton_in) const final {
        return SFBitsetToNHigherLevel(n_level, morton_in);
    }
    DefInt SFBitsetCoincideLevelVir(const DefSFBitset& morton_in) const final {
        return SFBitsetCoincideLevel(morton_in);
    }
    void SFBitsetFindAllNeighborsVir(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const final {
        std::array<DefSFBitset, 9> array_neighbors;
        SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
        ptr_vec_neighbors->resize(9);
        memcpy(ptr_vec_neighbors->data(), array_neighbors.data(),
            9 * sizeof(DefSFBitset));
    };


    std::vector<DefAmrLUint> GetMinBackgroundIndices() const final {
        return {k0MinIndexOfBackgroundNode_[kXIndex], k0MinIndexOfBackgroundNode_[kYIndex]};
    }
    std::vector<DefAmrLUint> GetMaxBackgroundIndices() const final {
        return {k0MaxIndexOfBackgroundNode_[kXIndex], k0MaxIndexOfBackgroundNode_[kYIndex]};
    }
    void GetMinM1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const final;
    void GetMaxP1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_p1_bitset) const final;
    void GetMinAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_bitsets) const final;
    void GetMaxAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_bitset) const final;
    void FindNeighboringCoarseFromFine(const DefSFBitset& sfbitset_fine,
        std::vector<DefSFBitset>* const ptr_vec_coarse) const final;

    SFBitsetAux2D() { SFBitsetInitial(); }

    // set protected members
    void SetSpaceBackground(const std::vector<DefReal>& space_background) override;

    // mpi related
 public:
    void GetNLevelCorrespondingOnes(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_last_ones) const override;
#ifdef ENABLE_MPI
    void FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefAmrLUint block_length,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min, const std::array<DefAmrLUint, 2>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefAmrLUint block_length,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min, const std::array<DefAmrLUint, 2>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefAmrLUint block_length,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min, const std::array<DefAmrLUint, 2>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefAmrLUint block_length,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min, const std::array<DefAmrLUint, 2>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;

 protected:
    void FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefAmrLUint block_length,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min, const std::array<DefAmrLUint, 2>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    bool CheckPartitionLimits(const DefSFBitset& code_in,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 2>& code_domain_min,
        const std::array<DefAmrLUint, 2>& code_domain_max) const;
#endif  // ENABLE_MPI
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
class  SFBitsetAux3D : public SFBitsetAuxInterface {
 public:
    // p is positive and n is negative
    static constexpr DefInt  kNodeIndexX0Y0Z0_ = 0, kNodeIndexXnY0Z0_ = 1,
     kNodeIndexXpY0Z0_ = 2, kNodeIndexX0YnZ0_ = 3, kNodeIndexX0YpZ0_ = 4,
     kNodeIndexX0Y0Zn_ = 9, kNodeIndexX0Y0Zp_ = 10, kNodeIndexXnYnZ0_ = 5,
     kNodeIndexXnYpZ0_ = 6, kNodeIndexXpYnZ0_ = 7, kNodeIndexXpYpZ0_ = 8,
     kNodeIndexXnY0Zn_ = 11, kNodeIndexXnY0Zp_ = 12, kNodeIndexXpY0Zn_ = 13,
     kNodeIndexXpY0Zp_ = 14, kNodeIndexX0YnZn_ = 15, kNodeIndexX0YnZp_ = 16,
     kNodeIndexX0YpZn_ = 17, kNodeIndexX0YpZp_ = 18, kNodeIndexXpYpZp_ = 19,
     kNodeIndexXpYnZp_ = 20, kNodeIndexXpYpZn_ = 21, kNodeIndexXpYnZn_ = 22,
     kNodeIndexXnYpZp_ = 23, kNodeIndexXnYnZp_ = 24, kNodeIndexXnYpZn_ = 25,
     kNodeIndexXnYnZn_ = 26;  ///< indices of node and  its 26 neighbors
    std::array<DefSFBitset, 3> SFBitsetMin_ = {~DefSFCodeToUint(0), ~DefSFCodeToUint(0), ~DefSFCodeToUint(0)},
        SFBitsetMax_ = {0, 0, 0};

    /* give an offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (k0MinIndexOfBackgroundNode_ * kDomainDx),
    and the default value is 1. */
    std::array<DefAmrLUint, 3> k0MinIndexOfBackgroundNode_ = {1, 1, 1};
    ///< the minimum index of background nodes in each direction*/
    std::array<DefAmrLUint, 3> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0MinIndexOfBackgroundNode_ is included in k0MaxIndexOfBackgroundNode_
    std::array<DefReal, 3> k0DomainSize_{};  ///< domain size
    std::array<DefReal, 3> k0DomainDx_{};  ///< grid space
    std::array<DefReal, 3> k0RealMin_{};  ///< k0MinIndexOfBackgroundNode_ * kDomainDx


    inline DefSFBitset SFBitsetToOneLowerLevel(const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToOneHigherLevel(const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNLowerLevel(const DefInt n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset SFBitsetToNHigherLevel(const DefInt n_level, const DefSFBitset& morton_in) const;
    inline DefSFBitset FindXNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindXPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindYPos(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindZNeg(const DefSFBitset& bitset_in) const;
    inline DefSFBitset FindZPos(const DefSFBitset& bitset_in) const;
    inline DefInt SFBitsetCoincideLevel(const DefSFBitset& bitset_in) const;

    DefSFBitset SFBitsetEncoding(const std::array<DefAmrLUint, 3>& coordi_index) const;
    void SFBitsetComputeIndices(const DefSFBitset& bitset_in, std::array<DefAmrLUint, 3>* const ptr_indices) const;
    void SFBitsetComputeCoordinate(const DefSFBitset& bitset_in,
    const std::array<DefReal, 3>& grid_space, std::array<DefReal, 3>* const ptr_coordi) const;

    DefSFBitset SFBitsetBitsForRefinement(const DefInt i_level) const;
    void SFBitsetSetMinAndMaxBounds(
       const std::array<DefAmrLUint, 3>& indices_min, const std::array<DefAmrLUint, 3>& indices_max);
    void SFBitsetNotOnDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::array<DefSFBitset, 3>& sfbitset_min,
        const std::array<DefSFBitset, 3>& sfbitset_max,
        std::array<bool, 3>* const ptr_bool_not_at_boundary_neg,
        std::array<bool, 3>* const ptr_bool_not_at_boundary_pos) const;
    void SFBitsetFindAllNeighbors(const DefSFBitset& sfbitset_center,
        std::array<DefSFBitset, 27>* const ptr_bitset_neighbour) const;
    void SFBitsetFindCellNeighbors(const DefSFBitset& sfbitset_corner,
        std::array<DefSFBitset, 8>* const ptr_vec_bitset_neighbour) const;
    template<typename DataType>
    bool SFBitsetBelongToOneCell(const DefSFBitset& sfbitset_in,
        const DefMap<DataType>& map_node_exist, std::array<DefSFBitset, 8>* const ptr_sfbitsets) const;
    bool SFBitsetBelongToOneCellAcrossTwoLevels(const DefSFBitset& sfbitset_in,
        const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_node_exist_coarse,
        std::array<std::pair<DefSFBitset, DefInt>, 8>* const ptr_sfbitsets) const;
    int ResetIndicesExceedingDomain(const std::array<DefAmrLUint, 3>& domain_min_indices,
        const std::array<DefAmrLUint, 3>& domain_max_indices,
        DefSFCodeToUint* const ptr_i_code, DefSFBitset* ptr_bitset_tmp) const;

    // virtual functions
    void SFBitsetInitial() override;
    DefSFBitset SFBitsetEncodingCoordi(
        const std::vector<DefReal>& grid_space, const std::vector<DefReal>& coordi) const override;
    void SFBitsetFindEdgeNode(const DefSFBitset& morton_in,
        std::vector<DefSFBitset>* const ptr_vec_edge_nodes) const final {
        ptr_vec_edge_nodes->resize(6);
        ptr_vec_edge_nodes->at(0) = FindXNeg(morton_in);
        ptr_vec_edge_nodes->at(1) = FindXPos(morton_in);
        ptr_vec_edge_nodes->at(2) = FindYNeg(morton_in);
        ptr_vec_edge_nodes->at(3) = FindYPos(morton_in);
        ptr_vec_edge_nodes->at(4) = FindZNeg(morton_in);
        ptr_vec_edge_nodes->at(5) = FindZPos(morton_in);
    }
    DefSFBitset FindNodeInNeg(const DefInt dir, const DefSFBitset& sfbitset_in) const override {
        if (dir == kXIndex) {
            return FindXNeg(sfbitset_in);
        } else if (dir == kYIndex) {
            return FindYNeg(sfbitset_in);
        } else if (dir == kZIndex) {
            return FindZNeg(sfbitset_in);
        } else {
             return sfbitset_in;
        }
    }
    DefSFBitset FindNodeInPos(const DefInt dir, const DefSFBitset& sfbitset_in) const override {
        if (dir == kXIndex) {
            return FindXPos(sfbitset_in);
        } else if (dir == kYIndex) {
            return FindYPos(sfbitset_in);
        } else if (dir == kZIndex) {
            return FindZPos(sfbitset_in);
        } else {
             return sfbitset_in;
        }
    }
    void SFBitsetComputeCoordinateVir(const DefSFBitset& bitset_in,
        const std::vector<DefReal>& grid_space, std::vector<DefReal>* const ptr_coordi) const final {
        std::array<DefReal, 3> arr_grid_space =
            {grid_space[kXIndex], grid_space[kYIndex], grid_space[kZIndex]}, arr_coordi;
        SFBitsetComputeCoordinate(bitset_in, arr_grid_space, &arr_coordi);
        ptr_coordi->resize(3);
        ptr_coordi->at(kXIndex) = arr_coordi[kXIndex];
        ptr_coordi->at(kYIndex) = arr_coordi[kYIndex];
        ptr_coordi->at(kZIndex) = arr_coordi[kZIndex];
    }
    void SFBitsetFindCellNeighborsVir(const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_bitset_neighbour) const final {
        std::array<DefSFBitset, 8> array_neighbors;
        SFBitsetFindCellNeighbors(sfbitset_corner, &array_neighbors);
        ptr_bitset_neighbour->resize(8);
        memcpy(ptr_bitset_neighbour->data(), array_neighbors.data(),
            8 * sizeof(DefSFBitset));
    }
    void FindNodesInReginOfGivenLength(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_nodes) const final;
    DefInt FindNodesInPeriodicRegionCorner(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    DefInt FindNodesInPeriodicRegionCornerOverlap(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const final;
    DefInt FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const DefInt region_length_neg, const DefInt region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    void FindNodesInPeriodicRegionCenter(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node) const final;
    DefInt FindNodesInPeriodicRegionCenterOverlap(const DefSFBitset& sfbitset_in,
        const std::vector<DefInt>& region_length_neg, const std::vector<DefInt>& region_length_pos,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_sfbitset_node,
        std::vector<std::pair<DefAmrLUint, DefSFBitset>>* const ptr_sfbitset_node_overlap) const final;
    DefSFBitset SetIindexinGivenDirection(const DefInt i_dir, const DefSFBitset& sfbitset_in,
        const DefSFBitset& sfbitset_target) const final;
    bool CheckExistenceCurrentLevel(
        const DefSFBitset& sfbitset_in, const DefMap<DefInt>& exist_nodes) const final;
    bool CheckIfOnGivenBoundary(const DefInt i_dir,
        const DefSFBitset& sfbitset_in, const DefSFBitset& sfbitset_boundary) const final;
    void FindNeighboringCoarseFromFine(const DefSFBitset& sfbitset_fine,
        std::vector<DefSFBitset>* const ptr_vec_coarse) const final;
    void SFBitsetHigherLevelInACell(const DefInt level_diff, const DefSFBitset& sfbitset_corner,
        std::vector<DefSFBitset>* const ptr_ptr_sfbitsets_higher_level) const final;
    void SFBitsetFindAllBondedNeighborsVir(const DefSFBitset& bitset_in,
        const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level,
        const std::vector<DefSFBitset>& domain_max_n_level,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const final;

    DefSFBitset  SFBitsetToNLowerLevelVir(const DefInt n_level, const DefSFBitset& morton_in) const final {
        return SFBitsetToNLowerLevel(n_level, morton_in);
    }
    DefSFBitset  SFBitsetToNHigherLevelVir(const DefInt n_level, const DefSFBitset& morton_in) const final {
        return SFBitsetToNHigherLevel(n_level, morton_in);
    }
    DefInt SFBitsetCoincideLevelVir(const DefSFBitset& morton_in) const final {
        return SFBitsetCoincideLevel(morton_in);
    }
    void SFBitsetFindAllNeighborsVir(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const final {
        std::array<DefSFBitset, 27> array_neighbors;
        SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
        ptr_vec_neighbors->resize(27);
        memcpy(ptr_vec_neighbors->data(), array_neighbors.data(),
            27 * sizeof(DefSFBitset));
    };

    std::vector<DefAmrLUint> GetMinBackgroundIndices() const final {
        return {k0MinIndexOfBackgroundNode_[kXIndex], k0MinIndexOfBackgroundNode_[kYIndex],
            k0MinIndexOfBackgroundNode_[kZIndex]};
    }
    std::vector<DefAmrLUint> GetMaxBackgroundIndices() const final {
        return {k0MaxIndexOfBackgroundNode_[kXIndex], k0MaxIndexOfBackgroundNode_[kYIndex],
            k0MaxIndexOfBackgroundNode_[kZIndex]};
    }
    void GetMinM1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const final;
    void GetMaxP1AtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_p1_bitset) const final;
    void GetMinAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_min,
        std::vector<DefSFBitset>* const ptr_min_bitsets) const final;
    void GetMaxAtGivenLevel(const DefInt i_level,
        std::vector<DefAmrLUint> indices_max,
        std::vector<DefSFBitset>* const ptr_max_bitset) const final;

    SFBitsetAux3D() { SFBitsetInitial(); }

    // set protected members
    void SetSpaceBackground(const std::vector<DefReal>& space_background) override;

    // mpi related
 public:
    void GetNLevelCorrespondingOnes(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_last_ones) const override;
#ifdef ENABLE_MPI
    void FindPartitionBlocksMax(const DefSFCodeToUint& code_in, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionBlocksMin(const DefSFCodeToUint& code_in, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionRemainMax(const DefSFCodeToUint& code_in, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartitionRemainMin(const DefSFCodeToUint& code_in, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;

 protected:
    void FindPartitionOneBlock(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition2BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition3BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition4BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition5BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition6BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition7BlocksMax(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition2BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition3BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition4BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition5BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition6BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    void FindPartition7BlocksMin(const DefSFBitset& bitset_corner, const DefInt block_level,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min, const std::array<DefAmrLUint, 3>& code_domain_max,
        DefMap<DefInt>* const ptr_map_partitioned) const;
    bool CheckPartitionLimits(const DefSFBitset& code_in,
        const DefSFCodeToUint& code_partition_min, const DefSFCodeToUint& code_partition_max,
        const std::array<DefAmrLUint, 3>& code_domain_min,
        const std::array<DefAmrLUint, 3>& code_domain_max) const;
#endif  // ENABLE_MPI
};
#endif  // DEBUG_DISABLE_3D_FUNCTION
}  // end namespace amrproject
}  // end namespace rootproject

// include inline implementation of morton code function
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
* @return  if true indicates the given node belongs to a cell
* @note ptr_sfbitsets[0]:(0, 0); ptr_sfbitsets[1]:(+x, 0);
*       ptr_sfbitsets[2]:(0, +y); ptr_sfbitsets[3]:(+x, +y);
*/
template<typename DataType>
bool SFBitsetAux2D::SFBitsetBelongToOneCell(
    const DefSFBitset& sfbitset_in, const DefMap<DataType>& map_node_exist,
    std::array<DefSFBitset, 4>* const ptr_sfbitsets)  const {
    if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
        return false;
    }
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1)) == map_node_exist.end()) {
        return false;
    }
    // (+x, +y)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    if (map_node_exist.find(ptr_sfbitsets->at(3)) == map_node_exist.end()) {
        return false;
    }
    // (0, +y)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2)) == map_node_exist.end()) {
        return false;
    }
    return true;
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to find space filling code constructing a cell (3D)
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  map_node_exist grid instance containing nodes at the same level
* @param[out] ptr_sfbitsets of nodes on the cell
* @return  if true indicates the given node belongs to a cell
* @note ptr_sfbitsets[0]:(0, 0, 0); ptr_sfbitsets[1]:(+x, 0, 0);
*       ptr_sfbitsets[2]:(0, +y, 0); ptr_sfbitsets[3]:(+x, +y, 0);
*       ptr_sfbitsets[4]:(0, 0, +z);  ptr_sfbitsets[5]:(+x, 0, +z);
*       ptr_sfbitsets[6]:(0, +y, +z); ptr_sfbitsets[7]: (+x, +y, +z).
*/
template<typename DataType>
bool SFBitsetAux3D::SFBitsetBelongToOneCell(
    const DefSFBitset& sfbitset_in, const DefMap<DataType>& map_node_exist,
    std::array<DefSFBitset, 8>* const ptr_sfbitsets) const {
    if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
        return false;
    }
    ptr_sfbitsets->at(0) = sfbitset_in;
    // (+x, 0, 0)
    ptr_sfbitsets->at(1) = FindXPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(1)) == map_node_exist.end()) {
        return false;
    }
    // (+x, +y, 0)
    ptr_sfbitsets->at(3) = FindYPos(ptr_sfbitsets->at(1));
    if (map_node_exist.find(ptr_sfbitsets->at(3))  == map_node_exist.end()) {
        return false;
    }
    // (0, +y, 0)
    ptr_sfbitsets->at(2) = FindYPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(2))  == map_node_exist.end()) {
        return false;
    }
    // (0, 0, +z)
    ptr_sfbitsets->at(4) = FindZPos(sfbitset_in);
    if (map_node_exist.find(ptr_sfbitsets->at(4))  == map_node_exist.end()) {
        return false;
    }
    // (+x, 0, +z)
    ptr_sfbitsets->at(5) = FindXPos(ptr_sfbitsets->at(4));
    if (map_node_exist.find(ptr_sfbitsets->at(5)) == map_node_exist.end()) {
        return false;
    }
    // (+x, +y, +z)
    ptr_sfbitsets->at(7) = FindYPos(ptr_sfbitsets->at(5));
    if (map_node_exist.find(ptr_sfbitsets->at(7)) == map_node_exist.end()) {
        return false;
    }
    // (0, +y, +z)
    ptr_sfbitsets->at(6) = FindYPos(ptr_sfbitsets->at(4));
    if (map_node_exist.find(ptr_sfbitsets->at(6)) == map_node_exist.end()) {
        return false;
    }
    return true;
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_GRID_SFBITSET_AUX_H_
