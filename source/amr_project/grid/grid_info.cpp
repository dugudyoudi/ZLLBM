//  Copyright (c) 2021 - 2023, Zhengliang Liu
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
}  // end namespace amrproject
}  // end namespace rootproject
