//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  define the class to manager other objects.
*/

#ifndef ROOTPROJECT_SOURCE_OBJ_MANAGER_H_
#define ROOTPROJECT_SOURCE_OBJ_MANAGER_H_
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
* @brief class used to manage the amr processess.
* @note using singleton design pattern
*/
class AmrManager {
 public:
    static AmrManager* GetInstance() {
        static AmrManager amr_instance_;
        return &amr_instance_;
    }

    // modules
#ifdef ENABLE_MPI
    std::shared_ptr<mpi::MpiManager> ptr_mpi_manager;
#endif
    std::shared_ptr<GridManagerInterface> ptr_grid_manager_;
    std::shared_ptr<IoManager> ptr_io_manager_;
    std::shared_ptr<CriterionManager> ptr_criterion_manager_;

    void LoadModules();

    void DefaultInitialization(DefUint dim, DefSizet level,
        int argc, char* argv[]);
    void SetupParameters();
    void InitialSimulation();
    void FinalizeSimulation();

    void SetTheSameLevelDependentInfoForAllLevels(
        std::shared_ptr<SolverCreatorInterface> ptr_solver_creator);

 private:
    static AmrManager* amr_instance_;
    AmrManager(void) {}
    ~AmrManager(void) {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_OBJ_MANAGER_H_