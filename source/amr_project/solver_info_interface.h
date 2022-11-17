//  Copyright (c) 2022, Zhengliang Liu
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
namespace grid {
class GridInfoCreatorInterface;
}  // end namespace grid
/**
* @class SolverInterface
* @brief abstract class used to manager solver.
*/
class SolverInterface {
public:
    DefUint k0SolverDims_ = 0;  ///< dimension
    std::string node_type_;
    std::shared_ptr<grid::GridInfoCreatorInterface>
        ptr_grid_info_creator_ = nullptr;
    virtual std::string GetSolverMethod() = 0;
    virtual void SolverInitial() = 0;
    virtual void RunSolver() = 0;
};
/**
* @class SolverCreatorInterface
* @brief abstract class used to generate sovler instance.
*/
class SolverCreatorInterface {
public:
    virtual std::shared_ptr<SolverInterface> CreateSolver() = 0;
};
}  // end namespace amrprproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
