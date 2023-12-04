//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file solver_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  define the class for Solver interface.
*/
#ifndef ROOTPROJECT_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
#define ROOTPROJECT_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
#include <string>
#include <memory>
namespace rootproject {
namespace amrproject {
class GridInfoCreatorInterface;
class GridInfoInterface;
class GridManagerInterface;
class SFBitsetAuxInterface;
/**
* @class SolverInterface
* @brief abstract class used to manager solver.
*/
class SolverInterface {
 public:
    DefAmrIndexUint k0SolverDims_ = 0;  ///< dimension
    std::string solver_type_;
    std::unique_ptr<GridInfoCreatorInterface> ptr_grid_info_creator_;
    GridManagerInterface* ptr_grid_manager_ = nullptr;
    virtual std::string GetSolverMethod() = 0;
    virtual void SetNodeFlagForSolver() = 0;
    virtual void SolverInitial() = 0;
    virtual void RunSolver(const DefReal time_step_current,
        const SFBitsetAuxInterface& sfbitset_aux,
        GridInfoInterface* const ptr_grid_info) = 0;
};
/**
* @class SolverCreatorInterface
* @brief abstract class used to generate solver instance.
*/
class SolverCreatorInterface {
 public:
    virtual std::shared_ptr<SolverInterface> CreateSolver() const = 0;
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
