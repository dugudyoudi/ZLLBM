//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/node_infor.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define parent class for nodes
*/

#ifndef ROOTPROJECT_SOURCE_LBM_LBM_INFOR_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_INFOR_H_
#include <string>
#include <vector>
#include "./defs_libs.h"
#include "./solver_info.h"
#include "grid/grid_info.h"
namespace rootproject {
namespace lbm {
/**
* @brief enumerate LBM collison operator type
*/
enum class ELbmCollisionOperatorType {
    kUndefined = 0,
    kLbmSrt =1,
    kLbmMrt =2
};
/**
* @brief  class to store LBM grid information at each level
*/
class GridInfoLBM final :public grid::GridInfoInterface {
 public:
    void initial_grid_node(const DefSFBitset& bit_set_in) override;
    void SetNumberOfVecElements() override;
};
class SolverLbmDnQn:public SolverInterface {
    friend class GridInfoLBM;
 public:
     DefAmrIndexUint k0NumQ_ = 0;
     ELbmCollisionOperatorType collision_type_ =
         ELbmCollisionOperatorType::kUndefined;
    std::string GetSolverMethod() final {
        return  "Solver LBM D" + std::to_string(k0SolverDims_)
            + "Q" + std::to_string(k0NumQ_) + " Model is active.";
    }
    void SolverInitial() override;
    void RunSolver()override {};
    void StreamAndCollid();

 protected:
    bool bool_vec_forces_ = false;
    // parameter to set indices for std::vector vec_real in grid::GridNode.
    // The default value is 2, indicating the 0 to (kNumQ - 1) elements are
    // distribution functions after collision (f^{\tilt}), while the kNumQ
    // to (2*kNumQ - 1) elements are distribution functions after stream (f).
    DefAmrIndexUint num_distribution_func_set_ = 2;
    ///< the mutiple of k0NumQ distribution functions need to be stored at each node.*/
    std::vector<DefReal> k0Cx_, k0Cy_, k0Cz_;
    ///< directions  of particle velocity
    std::vector<DefReal> k0Weights_;
    DefAmrIndexUint k0RhoIndex_ = 0, k0UxIndex_ = 0, k0UyIndex_ = 0, k0UzIndex_ = 0;
    DefAmrIndexUint k0FxIndex_ = 0, k0FyIndex_ = 0, k0FzIndex_ = 0;
};
class SolverLbmD2Q9 final :public SolverLbmDnQn {
    void SolverInitial() override;
};
class SolverCreatorLbmD2Q9 final :public SolverCreatorInterface {
    std::shared_ptr<SolverInterface> CreateSolver() override;
};
}  // end namespace lbm
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_INFOR_H_
