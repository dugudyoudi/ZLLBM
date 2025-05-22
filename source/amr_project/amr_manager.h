//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  define the class to manager other objects.
*/

#ifndef SOURCE_AMR_PROJECT_AMR_MANAGER_H_
#define SOURCE_AMR_PROJECT_AMR_MANAGER_H_
#include <memory>
#include <string>
#include <vector>
#include "../defs_libs.h"
#include "./solver_info_interface.h"
#include "io/io_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
#include "criterion/criterion_manager.h"
#include "io/input_parser.h"
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
    bool bool_write_checkpoint_ = false;

    // modules
#ifdef ENABLE_MPI
    int StartupAmrManager(int argc, char* argv[]);
#endif  // ENABLE_MPI
    std::unique_ptr<MpiManager> ptr_mpi_manager_;  // when MPI is not enabled, this point to an empty class

    std::unique_ptr<GridManagerInterface> ptr_grid_manager_;
    std::unique_ptr<IoManager> ptr_io_manager_;
    std::unique_ptr<CriterionManager> ptr_criterion_manager_;

    void LoadModules(DefInt dims);

    void StartupInitialization(DefInt dim);
    void SetupInputDependentParameters();
    void InitializeMesh();
    void InitializeAllSolvers();
    void RestartFromCheckPointMesh(DefInt time_step);
    void SetupMesh(DefInt time_step);
    void FinalizeSimulation();

    // time stepping related
    ETimeSteppingScheme k0TimeSteppingType_ = ETimeSteppingScheme::kMultiSteppingC2F;
    std::unique_ptr<TimeSteppingSchemeInterface> ptr_time_stepping_scheme_;
    void InstantiateTimeSteppingScheme();
    void TimeMarching(const DefAmrLUint time_step);

    void AddSolverToGridManager(const SolverCreatorInterface& solver_creator);
    void SetSameSolverDependentInfoForAllGrids(const std::shared_ptr<SolverInterface>& ptr_solver);

    void BroadCastInputParse(InputParser* const ptr_input_parser) const;

 private:
    DefAmrLUint current_time_step_ = 0;

    // mpi
    MpiManager* GetPointerToMpiManager() const {
        if (ptr_mpi_manager_ != nullptr) {
            return ptr_mpi_manager_.get();
        } else {
            return nullptr;
        }
    }

    static AmrManager* amr_instance_;
    AmrManager(void) {}
    ~AmrManager(void) {}

#ifdef DEBUG_CHECK_GRID
    // functions for the purpose of debug
 public:
    void CheckMeshAfterInitialization() const;
#endif
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_AMR_MANAGER_H_
