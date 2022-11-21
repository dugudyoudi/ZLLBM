//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file io.h
* @author Zhengliang Liu
* @date  2022-5-18
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
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
    static std::string logfile_name;

    void DefaultInitialization();
    void SetIoParameters();

    // write flow field options
    bool bool_binary_ = true;
    OutputDataFormat k0OutputDataFormat_;

    void OutputFlowfield(
        std::shared_ptr<GridManagerInterface> ptr_grid_manager,
        std::shared_ptr<CriterionManager> ptr_criterion_manager);
    //VtkWriterManager vtk_instance_;
};
}  // end amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_IO_IO_MANAGER_H_
