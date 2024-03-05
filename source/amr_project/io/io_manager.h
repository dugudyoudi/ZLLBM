//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file io_manager.h
* @author Zhengliang Liu
* @date  2022-5-18
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
#include <string>
#include "../defs_libs.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace amrproject {
/**
* @class IoManager
* @brief class used to manage io processes.
* @date  2022-5-20
*/
class IoManager {
 public:

    void DefaultInitialization();
    void SetIoParameters();

    // write flow field options
    bool bool_binary_ = true;
    OutputDataFormat k0OutputDataFormat_;

    void OutputFlowfield(const std::string& prog_name,
        GridManagerInterface* const ptr_grid_manager,
        CriterionManager* const ptr_criterion_manager);
    VtkWriterManager vtk_instance_;
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
