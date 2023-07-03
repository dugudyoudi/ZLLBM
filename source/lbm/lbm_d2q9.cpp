//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include "lbm/lbm_info.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbm {
std::shared_ptr<SolverInterface> SolverCreatorLbmD2Q9::CreateSolver() {
    std::shared_ptr<SolverLbmD2Q9> ptr_temp =
        std::make_shared<lbm::SolverLbmD2Q9>();
    ptr_temp->k0SolverDims_ = 2;
    ptr_temp->k0NumQ_ = 9;
    ptr_temp->node_type_ = "LbmD2Q9";
    return ptr_temp;
}
}  // end namespace lbm
}  // end namespace rootproject