//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
*/
#include <string>
#include <mpi.h>
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
    int flag_node = 0;
    // x min
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMin_[kXIndex].to_ullong();
    if (code_one_dim <= code_boundary) {
        if (code_one_dim == code_boundary) {
            flag_node |= 1;
        } else {
            return -1;
        }
    }
    // x max
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeXRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMax_[kXIndex].to_ullong();
    if (code_one_dim >= code_boundary) {
        if (code_one_dim == code_boundary) {
            flag_node |= 8;
        } else {
            return -1;
        }
    }
    // y min
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMin_[kYIndex].to_ullong();
    if (code_one_dim <= code_boundary) {
        if (code_one_dim == code_boundary) {
            flag_node |= 2;
        } else {
            return -1;
        }
    }
    // y max
    code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeYRef_[current_bit]).to_ullong();
    code_boundary = k0VecBitsetDomainMax_[kYIndex].to_ullong();
    if (code_one_dim >= code_boundary) {
        if (code_one_dim == code_boundary) {
            flag_node |= 16;
        } else {
            return -1;
        }
    }
    if (dims == 3) {
        // z min
        code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]).to_ullong();
        code_boundary = k0VecBitsetDomainMin_[kZIndex].to_ullong();
        if (code_one_dim <= code_boundary) {
            if (code_one_dim == code_boundary) {
                flag_node |= 4;
            } else {
                return -1;
            }
        }
        // z max
        code_one_dim = (bitset_in&sfbitset_aux.k0SFBitsetTakeZRef_[current_bit]).to_ullong();
        code_boundary = k0VecBitsetDomainMax_[kZIndex].to_ullong();
        if (code_one_dim >= code_boundary) {
            if (code_one_dim == code_boundary) {
                flag_node |= 32;
            } else {
                return -1;
            }
        }
    }
    return flag_node;
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
