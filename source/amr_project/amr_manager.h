//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  define the class to manager other objects.
*/

#ifndef ROOTPROJECT_SOURCE_OBJ_MANAGER_H_
#define ROOTPROJECT_SOURCE_OBJ_MANAGER_H_
#include <memory>
#include <string>
#include <vector>
#include "../defs_libs.h"
#include "./solver_info_interface.h"
#include "io/io_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @class AmrManager
* @brief class used to manage the amr processes.
* @note using singleton design pattern
*/
class AmrManager {
 public:
    std::string program_name_;
    static AmrManager* GetInstance() {
        static AmrManager amr_instance_;
        return &amr_instance_;
    }

    // modules
#ifdef ENABLE_MPI
    int SetUpProgramFeature(int argc, char* argv[]);
    std::unique_ptr<MpiManager> ptr_mpi_manager_;
#endif  // ENABLE_MPI
    std::unique_ptr<GridManagerInterface> ptr_grid_manager_;
    std::unique_ptr<IoManager> ptr_io_manager_;
    std::unique_ptr<CriterionManager> ptr_criterion_manager_;

    void LoadModules(DefAmrIndexUint dims);

    void DefaultInitialization(DefAmrIndexUint dim, DefAmrIndexUint level);
    void SetupParameters();
    void InitializeMesh();
    void SetupSolverForGrids();
    void FinalizeSimulation();

    // time stepping related
    ETimeSteppingScheme k0TimeSteppingType_ = ETimeSteppingScheme::kMultiSteppingC2F;
    std::unique_ptr<TimeSteppingSchemeInterface> ptr_time_stepping_scheme_;
    void InstantiateTimeSteppingScheme();
    void TimeMarching(const DefAmrIndexLUint time_step);

    void AddSolverToGridManager(const SolverCreatorInterface& solver_creator);
    void SetDependentInfoForAllLevelsTheSame(const std::shared_ptr<SolverInterface>& ptr_solver);

 private:
    static AmrManager* amr_instance_;
    AmrManager(void) {}
    ~AmrManager(void) {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_OBJ_MANAGER_H_