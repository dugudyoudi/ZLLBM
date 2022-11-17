//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file io_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage io processes.
* @date  2022-8-6
*/
#include "io/log_write.h"
#include "io/io_manager.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject{
namespace io {
/**
* @brief function to setup default io related parameters.
*/
void IoManager::DefaultInitialization() {
    // set header of the log file (start time)
    logfile_name = "log_of_node";
    LogStartTime();

#ifdef ENABLE_MPI
    int rank_id = 0;
    rank_id = mpi::MpiManager::rank_id_;
    if (rank_id == 0) io::LogInfo(
        "number of processes is: "
        + std::to_string(mpi::MpiManager::GetInstance()->numb_of_ranks_));
#endif  // ENABLE_MPI
    bool_binary_ = true;

    // set format for writing output data
    if (typeid(k0OutputDataFormat_.output_real_.cast_type(1.))
        == typeid(float)) {
        k0OutputDataFormat_.output_real_.printf_format_ = "%f";
        k0OutputDataFormat_.output_real_.format_name_ = "Float32";
    } else if (typeid(k0OutputDataFormat_.output_real_.cast_type(1.))
        == typeid(double)) {
        k0OutputDataFormat_.output_real_.printf_format_ = "%lf";
        k0OutputDataFormat_.output_real_.format_name_ = "Float64";
    } else {
        LogWarning("Output format for real type was not found.");
    }

    if (typeid(k0OutputDataFormat_.output_uint_.cast_type(1))
        == typeid(uint16_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%hu";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt16";
    } else if (typeid(k0OutputDataFormat_.output_uint_.cast_type(1))
        == typeid(uint32_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%u";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt32";
    } else if (typeid(k0OutputDataFormat_.output_uint_.cast_type(1))
        == typeid(uint64_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%llu";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt64";
    } else {
        LogWarning("Output format for uint type was not found.");
    }

    if (typeid(k0OutputDataFormat_.output_int_.cast_type(1))
        == typeid(int16_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%hd";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt16";
    } else if (typeid(k0OutputDataFormat_.output_uint_.cast_type(1))
        == typeid(int32_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%d";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt32";
    } else if (typeid(k0OutputDataFormat_.output_uint_.cast_type(1))
        == typeid(int64_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%lld";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt64";
    } else {
        LogWarning("Output format for int type was not found.");
    }
}
/**
* @brief function to set io related parameters.
*/
void IoManager::SetIoParameters() {
    //vtk_instance_.OptionInitial(bool_binary_);
}
/**
* @brief function to write flowfield.
*/
void IoManager::OutputFlowfield(
    std::shared_ptr<grid::GridManagerInterface> ptr_grid_manager,
    std::shared_ptr<criterion::CriterionManager> ptr_criterion_manager){
    //vtk_instance_.WriteVtu(bool_binary_, k0OutputDataFormat_,
    //    ptr_grid_manager, ptr_criterion_manager);
}
}  // end namespace io
}  // end amrproject
}  // end namespace rootproject