//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_interface.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
#include <string>
#include <vector>
#include <memory>
#include "../defs_libs.h"
#include "solver_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace lbmproject {
class SolverCreatorLbmD2Q9;
/**
* @brief enumerate LBM collision operator type
*/
enum class ELbmCollisionOperatorType {
    kUndefined = 0,
    kLbmSrt = 1,
    kLbmSrtForce = 2,
    kLbmMrt = 3,
    kLbmMrtForce = 4
};
struct GridNodeLbm : public amrproject::GridNode {
    std::vector<DefReal> f_collide_, f_;
    DefReal rho_;
    std::vector<DefReal> velocity_;
    std::vector<DefReal> force_;
};
/**
* @brief  class to store LBM grid information at each refinement level
*/
class GridInfoLbmInteface : public amrproject::GridInfoInterface {
 public:
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_ = nullptr;
    void SetPointerToCurrentNodeType() override;
    DefMap<std::unique_ptr<GridNodeLbm>>* GetPointerToLbmGrid();
    std::unique_ptr<amrproject::GridNode> GridNodeCreator() override {
        return std::make_unique<GridNodeLbm>();
    }
    void InitialGridNode(const DefSFBitset& bitset_in) override;
};
class GridInfoLbmIntefaceCreator :
    public amrproject::GridInfoCreatorInterface {
 public:
    std::shared_ptr<amrproject::GridInfoInterface>
        CreateGridInfo() override {
        return std::make_shared<GridInfoLbmInteface>();
    };
};
class SolverLbmInterface :public amrproject::SolverInterface {
    friend class SolverCreatorLbmD2Q9;

 public:
    DefAmrIndexUint k0NumQ_ = 0;
    ELbmCollisionOperatorType collision_type_ = ELbmCollisionOperatorType::kUndefined;
    std::string GetSolverMethod() final {
        return  "Solver LBM D" + std::to_string(k0SolverDims_)
            + "Q" + std::to_string(k0NumQ_) + " Model is active.";
    }

    SolverLbmInterface() {
        ptr_grid_info_creator_ = std::make_unique<GridInfoLbmIntefaceCreator>();
    }
    void SolverInitial() final;
    void RunSolver() override {};

    virtual void InitialSetIndices() = 0;
    virtual void Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAux2D& sfbitset_aux2d,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const = 0;
    void StreamAndCollide() const;

 protected:
    bool bool_forces_ = false;
    std::vector<DefReal> k0Cx_, k0Cy_, k0Cz_;
    ///< directions  of particle velocity
    std::vector<DefReal> k0Weights_;
    void (SolverLbmInterface::*ptr_func_cal_feq_)(
        const GridNodeLbm&, std::vector<DefReal>* const) const = nullptr;
    inline DefReal Square(DefReal x) const {return x * x;}
    ///< the mutiple of k0NumQ distribution functions need to be stored at each node.*/
    void CalFeq2D(const GridNodeLbm& node, std::vector<DefReal>* const ptr_feq) const;
    void CalFeq3D(const GridNodeLbm& node, std::vector<DefReal>* const ptr_feq) const;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
