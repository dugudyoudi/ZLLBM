//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/node_infor.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define parent class for nodes
*/

#ifndef ROOTPROJECT_SOURCE_LBM_LBM_MANAGER_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_MANAGER_H_
#include <memory>
#include <vector>
#include "./defs_libs.h"
#include "lbm/lbm_info.h"
#include "grid/grid_info.h"
namespace rootproject {
namespace lbm {
class LbmManager {
public:
   DefUint  k0LbmDims_;  ///< dimension
};
}  // end namespace lbm
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_MANAGER_H_
