//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file solver_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  define the class for Solver interface.
*/
#ifndef SOURCE_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
#define SOURCE_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
#include <string>
#include <memory>
#include <vector>
#include "grid/grid_enumerates.h"
#include "io/input_parser.h"
namespace rootproject {
namespace amrproject {
class GridInfoCreatorInterface;
class GridInfoInterface;
class GridManagerInterface;
class SFBitsetAuxInterface;
class MpiManager;
/**
* @class SolverInterface
* @brief abstract class used to manager solver.
*/
class SolverInterface {
 protected:
    DefInt k0SolverDims_ = 0;  ///< dimension
    std::string solver_type_, name_;  // name_ must be unique for each solver instance when using inputparser
    GridManagerInterface* ptr_grid_manager_ = nullptr;

 public:
    // set and get functions
    void SetSolverDims(const DefInt k0SolverDims) {k0SolverDims_ = k0SolverDims;}
    void SetName(const std::string& name) {name_ = name;}
    void SetPtrToGridManager(GridManagerInterface* const ptr_grid_manager) {
        ptr_grid_manager_ = ptr_grid_manager;
    }
    DefInt GetSolverDim() const {return k0SolverDims_;}
    std::string GetSolverType() const {return solver_type_;}
    std::string GetName() const {return name_;}
    GridManagerInterface* GetPtrToParentGridManager() const {return ptr_grid_manager_;}

    std::unique_ptr<GridInfoCreatorInterface> ptr_grid_info_creator_;
    virtual std::string GetSolverMethod() = 0;
    virtual void SolverInitial() = 0;
    virtual void RunSolverOnGivenGrid(const ETimeSteppingScheme time_scheme,
        const DefInt time_step_level, const DefReal time_step_current,
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) = 0;
    virtual void InformationFromGridOfDifferentLevel(
        const DefInt time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) {}
    virtual void ReadAndSetupSolverParameters(const InputParser& input_parser) {}

    virtual ~SolverInterface() = default;
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
#endif  // SOURCE_AMR_PROJECT_SOLVER_INFO_INTERFACE_H_
