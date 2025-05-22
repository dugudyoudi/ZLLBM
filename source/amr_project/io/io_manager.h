//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file io_manager.h
* @author Zhengliang Liu
* @brief define the class used to manage IO processes for grids.
* @date  2022-5-18
*/
#ifndef SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
#define SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
#include <string>
#include "../defs_libs.h"
#include "io/vtk_writer.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @class IoManager
* @brief class used to manage io processes.
* @date  2022-5-20
*/
class IoManager {
 public:
    void SetupOutputFormat();
    void SetupDependentIOParameters();

    // write output data options
    bool bool_binary_ = true;
    OutputDataFormat k0OutputDataFormat_;

    void OutputMeshData(const std::string& prog_name,
        GridManagerInterface* const ptr_grid_manager,
        CriterionManager* const ptr_criterion_manager);
    VtkWriterManager vtk_instance_;

    // write checkpoint data
    void WriteCheckPointData(const DefAmrLUint time_step, const std::string& prog_name,
        const SFBitsetAuxInterface& sfbitset_aux,
        const GridManagerInterface& grid_manager,
        const CriterionManager& criterion_manager, const MpiManager& mpi_manager) const;
    void ReadCheckPointData(const DefAmrLUint time_step, const std::string& prog_name,
        GridManagerInterface* const ptr_grid_manager,
        CriterionManager* const ptr_criterion_manager, MpiManager* const mpi_manager) const;
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
