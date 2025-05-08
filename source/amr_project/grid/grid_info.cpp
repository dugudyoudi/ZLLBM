//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_info.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid information.
* @date  2022-6-7
*/
#include <string>
#include "grid/grid_info_interface.h"
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {

void GridInfoInterface::SetGridLevel(DefInt i_level) {
    i_level_ = i_level;
}
void GridInfoInterface::SetComputationalCost(DefInt computational_cost) {
    computational_cost_ = computational_cost;
}
void GridInfoInterface::SetGridSpace(const std::vector<DefReal>& grid_space) {
    if (static_cast<DefInt>(grid_space.size()) != GetPtrToParentGridManager()->k0GridDims_) {
        LogManager::LogError("Dimension of grid space is not the same as the dimension of grid.");
    } else {
        grid_space_ = grid_space;
    }
}
void GridInfoInterface::SetPtrSFBitsetAux(SFBitsetAuxInterface* const ptr_sfbitset_aux) {
    ptr_sfbitset_aux_ = ptr_sfbitset_aux;
}
void GridInfoInterface::SetPtrSolver(const std::weak_ptr<SolverInterface>& ptr_solver) {
    ptr_solver_ = ptr_solver;
}
void GridInfoInterface::SetNumFine2CoarseLayer(DefInt num_fine2coarse_layer) {
    k0NumFine2CoarseLayer_ = num_fine2coarse_layer;
}
void GridInfoInterface::SetNumCoarse2FineLayer(DefInt num_coarse2fine_layer) {
    k0NumCoarse2FineLayer_ = num_coarse2fine_layer;
}
/**
 * @brief function to set the value of a member variable identified by a string.
 * @param str_var the name of the member variable to be set
 * @param value the value to set the member variable to
 */
void GridInfoInterface::SetMemberVariable(const std::string& str_var, int value) {
    if (kMemberNames_.find(str_var) != kMemberNames_.end()) {
        void* ptr_member = kMemberNames_[str_var];
        int* ptr_typed_member = static_cast<int*>(ptr_member);

        if (typeid(*ptr_typed_member) == typeid(value)) {
            *ptr_typed_member = value;
        } else {
            LogManager::LogError("type mismatch: member variable " + str_var
                + " requires a value of type " + typeid(*ptr_typed_member).name());
        }
    } else {
        LogManager::LogError("variable " + str_var + " not found.");
    }
}
/**
 * @brief function to get computational domain infomation at current level.
 * @return domain information.
 */
DomainInfo GridInfoInterface::GetDomainInfo() const {
    DomainInfo domain_info;
    const std::vector<DefAmrLUint> indices_min = ptr_parent_grid_manager_->GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = ptr_parent_grid_manager_->GetMaxIndexOfBackgroundNodeArrAsVec();
    ptr_sfbitset_aux_->GetMinAtGivenLevel(i_level_, indices_min, &domain_info.domain_min_n_level_);
    ptr_sfbitset_aux_->GetMaxAtGivenLevel(i_level_, indices_max, &domain_info.domain_max_n_level_);
    domain_info.bool_periodic_domain_ = CheckIfPeriodicDomainRequired(ptr_parent_grid_manager_->k0GridDims_,
        &domain_info.periodic_min_, &domain_info.periodic_max_);
    domain_info.grid_space_ = GetGridSpace();
    return domain_info;
}
/**
 * @brief function to check if a node is inside, on or outside the cubic domain boundary.
 * @param[in] dims dimension of the grid.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @return flag indicate node on which domain boundaries.
 */
DefInt GridInfoInterface::CheckIfNodeOutsideCubicDomain(const DefInt dims,
    const DefSFBitset& bitset_in, const SFBitsetAuxInterface& sfbitset_aux) const {
    DefInt current_bit = sfbitset_aux.kRefCurrent_;
    DefSFCodeToUint code_one_dim, code_boundary;
    DefInt flag_node = kFlagInsideDomain_;
    const std::array<DefSFBitset, 2>& take_xref = sfbitset_aux.GetTakeXRef(),
        take_yref =  sfbitset_aux.GetTakeYRef(), take_zref = sfbitset_aux.GetTakeZRef();
    // x min
    code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_xref[current_bit]);
    code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMin_[kXIndex]);
    if (code_one_dim < code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagXMinBoundary_;
    }
    // x max
    code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_xref[current_bit]);
    code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMax_[kXIndex]);
    if (code_one_dim > code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagXMaxBoundary_;
    }
    // y min
    code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_yref[current_bit]);
    code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMin_[kYIndex]);
    if (code_one_dim < code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagYMinBoundary_;
    }
    // y max
    code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_yref[current_bit]);
    code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMax_[kYIndex]);
    if (code_one_dim > code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagYMaxBoundary_;
    }
    if (dims == 3) {
        // z min
        code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_zref[current_bit]);
        code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMin_[kZIndex]);
        if (code_one_dim < code_boundary) {
            return kFlagOutsideDomain_;
        } else if (code_one_dim == code_boundary) {
            flag_node |= kFlagZMinBoundary_;
        }
        // z max
        code_one_dim = sfbitset_aux.SFBitsetToSFCode(bitset_in&take_zref[current_bit]);
        code_boundary = sfbitset_aux.SFBitsetToSFCode(k0VecBitsetDomainMax_[kZIndex]);
        if (code_one_dim > code_boundary) {
            return kFlagOutsideDomain_;
        } else if (code_one_dim == code_boundary) {
            flag_node |= kFlagZMaxBoundary_;
        }
    }
    return flag_node;
}
/**
 * @brief function to check if there are nodes on periodic boundaries of a cubic domain.
 * @param[in] dims dimension of the grid.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[in] periodic_min booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_nodes_periodic pointer to counterpart nodes on the periodic boundaries.
 */
void GridInfoInterface::CheckNodesOnCubicPeriodicBoundary(const DefInt dims, const DefSFBitset& bitset_in,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const SFBitsetAuxInterface& sfbitset_aux, std::vector<DefSFBitset>* const ptr_nodes_periodic) const {
    ptr_nodes_periodic->clear();
    DefInt current_bit = sfbitset_aux.kRefCurrent_;
    DefInt other_bit = sfbitset_aux.kRefOthers_;
    const std::array<DefSFBitset, 2>& take_xref = sfbitset_aux.GetTakeXRef(),
        take_yref =  sfbitset_aux.GetTakeYRef(), take_zref = sfbitset_aux.GetTakeZRef();
    // x min
    if (periodic_min.at(kXIndex)
        && ((bitset_in&take_xref[current_bit]) == k0VecBitsetDomainMin_[kXIndex])) {
        ptr_nodes_periodic->emplace_back((bitset_in&take_xref[other_bit])
            |k0VecBitsetDomainMax_[kXIndex]);
    }
    // x max
    if (periodic_max.at(kXIndex)
        && ((bitset_in&take_xref[current_bit]) == k0VecBitsetDomainMax_[kXIndex])) {
        ptr_nodes_periodic->emplace_back((bitset_in&take_xref[other_bit])
            |k0VecBitsetDomainMin_[kXIndex]);
    }
    // y min
    if (periodic_min.at(kYIndex)
        && ((bitset_in&take_yref[current_bit]) == k0VecBitsetDomainMin_[kYIndex])) {
        ptr_nodes_periodic->emplace_back((bitset_in&take_yref[other_bit])
            |k0VecBitsetDomainMax_[kYIndex]);
    }
    // y max
    if (periodic_max.at(kYIndex)
        && ((bitset_in&take_yref[current_bit]) == k0VecBitsetDomainMax_[kYIndex])) {
        ptr_nodes_periodic->emplace_back((bitset_in&take_yref[other_bit])
            |k0VecBitsetDomainMin_[kYIndex]);
    }
    if (dims == 3) {
        // z min
        if (periodic_min.at(kZIndex)
            && ((bitset_in&take_zref[current_bit]) == k0VecBitsetDomainMin_[kZIndex])) {
            ptr_nodes_periodic->emplace_back((bitset_in&take_zref[other_bit])
                |k0VecBitsetDomainMax_[kZIndex]);
        }
        // z max
        if (periodic_max.at(kZIndex)
            && ((bitset_in&take_zref[current_bit]) == k0VecBitsetDomainMax_[kZIndex])) {
            ptr_nodes_periodic->emplace_back((bitset_in&take_zref[other_bit])
                |k0VecBitsetDomainMin_[kZIndex]);
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
