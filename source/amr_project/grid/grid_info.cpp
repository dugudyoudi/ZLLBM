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
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
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
 * @brief function to check if a node is inside, on or outside the cubic domain boundary.
 * @param[in] dims dimension of the grid.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @return flag indicate node on which domain boundaries, 1: x min, 8: x max, 2: y min, 16: y max, 4: z min, 32: z max,
 *         0 inside domain, -1 outside domain
 */
int GridInfoInterface::CheckIfNodeOutsideCubicDomain(const DefAmrIndexUint dims,
    const DefSFBitset& bitset_in, const SFBitsetAuxInterface& sfbitset_aux) const {
    DefAmrIndexUint current_bit = sfbitset_aux.kRefCurrent_;
    DefSFCodeToUint code_one_dim, code_boundary;
    int flag_node = kFlagInsideDomain_;
    // x min
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMin_[kXIndex].to_ullong();
    if (code_one_dim < code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagXMinBoundary_;
    }
    // x max
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMax_[kXIndex].to_ullong();
    if (code_one_dim > code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagXMaxBoundary_;
    }
    // y min
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMin_[kYIndex].to_ullong();
    if (code_one_dim < code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagYMinBoundary_;
    }
    // y max
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMax_[kYIndex].to_ullong();
    if (code_one_dim > code_boundary) {
        return kFlagOutsideDomain_;
    } else if (code_one_dim == code_boundary) {
        flag_node |= kFlagYMaxBoundary_;
    }
    if (dims == 3) {
        // z min
        code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]).to_ullong();
        code_boundary = k0VecBitsetDomainMin_[kZIndex].to_ullong();
        if (code_one_dim < code_boundary) {
            return kFlagOutsideDomain_;
        } else if (code_one_dim == code_boundary) {
            flag_node |= kFlagZMinBoundary_;
        }
        // z max
        code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]).to_ullong();
        code_boundary = k0VecBitsetDomainMax_[kZIndex].to_ullong();
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
void GridInfoInterface::CheckNodesOnCubicPeriodicBoundary(const DefAmrIndexUint dims, const DefSFBitset& bitset_in,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const SFBitsetAuxInterface& sfbitset_aux, std::vector<DefSFBitset>* const ptr_nodes_periodic) const {
    ptr_nodes_periodic->clear();
    DefAmrIndexUint current_bit = sfbitset_aux.kRefCurrent_;
    DefAmrIndexUint other_bit = sfbitset_aux.kRefOthers_;
    // x min
    if (periodic_min.at(kXIndex)
        && ((bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]) == k0VecBitsetDomainMin_[kXIndex])) {
        ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[other_bit])
            |k0VecBitsetDomainMax_[kXIndex]);
    }
    // x max
    if (periodic_max.at(kXIndex)
        && ((bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]) == k0VecBitsetDomainMax_[kXIndex])) {
        ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[other_bit])
            |k0VecBitsetDomainMin_[kXIndex]);
    }
    // y min
    if (periodic_min.at(kYIndex)
        && ((bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]) == k0VecBitsetDomainMin_[kYIndex])) {
        ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[other_bit])
            |k0VecBitsetDomainMax_[kYIndex]);
    }
    // y max
    if (periodic_max.at(kYIndex)
        && ((bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]) == k0VecBitsetDomainMax_[kYIndex])) {
        ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[other_bit])
            |k0VecBitsetDomainMin_[kYIndex]);
    }
    if (dims == 3) {
        // z min
        if (periodic_min.at(kZIndex)
            && ((bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]) == k0VecBitsetDomainMin_[kZIndex])) {
            ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[other_bit])
                |k0VecBitsetDomainMax_[kZIndex]);
        }
        // z max
        if (periodic_max.at(kZIndex)
            && ((bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]) == k0VecBitsetDomainMax_[kZIndex])) {
            ptr_nodes_periodic->push_back((bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[other_bit])
                |k0VecBitsetDomainMin_[kZIndex]);
        }
    }
}
/**
 * @brief function to check if information from coarse node is needed.
 * @param timing enum class to identify timing in current time step.
 * @param time_scheme enum class to identify time stepping scheme used in computation.
 * @param time_step count of time steps of current level.
 * @return true: need for information from coarse grid.
 */
bool GridInfoInterface::CheckNeedInfoFromCoarse(const ETimingInOneStep timing,
    const ETimeSteppingScheme time_scheme, const DefAmrIndexUint time_step) const {
    switch (time_scheme) {
    case ETimeSteppingScheme::kMultiSteppingC2F: {
        if (timing == ETimingInOneStep::kStepEnd) {
            if (time_step%2 == 0) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
        break;
    }
    default:
        LogManager::LogError("undefined time stepping type in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        break;
    }
    return false;
}
/**
 * @brief function to check if information from fine node is needed.
 * @param timing enum class to identify timing in current time step.
 * @param time_scheme enum class to identify time stepping scheme used in computation.
 * @param time_step count of time steps of current level.
 * @return true: need for information from fine grid.
 */
bool GridInfoInterface::CheckNeedInfoFromFine(const ETimingInOneStep timing,
    const ETimeSteppingScheme time_scheme, const DefAmrIndexUint time_step) const {
    switch (time_scheme) {
    case ETimeSteppingScheme::kMultiSteppingC2F: {
        if (timing == ETimingInOneStep::kStepBegin) {
            return true;
        } else {
            return false;
        }
        break;
    }
    default:
        LogManager::LogError("undefined time stepping type in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        break;
    }
    return false;
}
}  // end namespace amrproject
}  // end namespace rootproject
