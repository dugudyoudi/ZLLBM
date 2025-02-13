//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file io_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage io processes.
* @date  2022-8-6
*/
#include <string>
#include "io/log_write.h"
#include "io/io_manager.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to setup default io related parameters.
*/
void IoManager::SetupOutputFormat() {
    // set format for writing output data
    if (typeid(k0OutputDataFormat_.output_real_.CastType(1.))
        == typeid(float)) {
        k0OutputDataFormat_.output_real_.printf_format_ = "%f";
        k0OutputDataFormat_.output_real_.format_name_ = "Float32";
    } else if (typeid(k0OutputDataFormat_.output_real_.CastType(1.))
        == typeid(double)) {
        k0OutputDataFormat_.output_real_.printf_format_ = "%lf";
        k0OutputDataFormat_.output_real_.format_name_ = "Float64";
    } else {
        LogManager::LogWarning("Output format for real type was not found in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    if (typeid(k0OutputDataFormat_.output_uint_.CastType(1))
        == typeid(uint16_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%hu";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt16";
    } else if (typeid(k0OutputDataFormat_.output_uint_.CastType(1))
        == typeid(uint32_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%u";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt32";
    } else if (typeid(k0OutputDataFormat_.output_uint_.CastType(1))
        == typeid(uint64_t)) {
        k0OutputDataFormat_.output_uint_.printf_format_ = "%llu";
        k0OutputDataFormat_.output_uint_.format_name_ = "UInt64";
    } else {
        LogManager::LogWarning("Output format for uint type was not found in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    if (typeid(k0OutputDataFormat_.output_sizet_.CastType(1))
        == typeid(uint16_t)) {
        k0OutputDataFormat_.output_sizet_.printf_format_ = "%hu";
        k0OutputDataFormat_.output_sizet_.format_name_ = "UInt16";
    } else if (typeid(k0OutputDataFormat_.output_sizet_.CastType(1))
        == typeid(uint32_t)) {
        k0OutputDataFormat_.output_sizet_.printf_format_ = "%u";
        k0OutputDataFormat_.output_sizet_.format_name_ = "UInt32";
    } else if (typeid(k0OutputDataFormat_.output_sizet_.CastType(1))
        == typeid(uint64_t)) {
        k0OutputDataFormat_.output_sizet_.printf_format_ = "%llu";
        k0OutputDataFormat_.output_sizet_.format_name_ = "UInt64";
    } else {
        LogManager::LogWarning("Output format for uint type was not found in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    if (typeid(k0OutputDataFormat_.output_int_.CastType(1))
        == typeid(int16_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%hd";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt16";
    } else if (typeid(k0OutputDataFormat_.output_uint_.CastType(1))
        == typeid(int32_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%d";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt32";
    } else if (typeid(k0OutputDataFormat_.output_uint_.CastType(1))
        == typeid(int64_t)) {
        k0OutputDataFormat_.output_int_.printf_format_ = "%lld";
        k0OutputDataFormat_.output_int_.format_name_ = "UInt64";
    } else {
        LogManager::LogWarning("Output format for int type was not found in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief function to set io related parameters.
*/
void IoManager::SetupDependentIOParameters() {
    vtk_instance_.OptionInitial(bool_binary_);
}
/**
* @brief function to write flow field.
*/
void IoManager::OutputFlowField(
    const std::string& prog_name,
    GridManagerInterface* const ptr_grid_manager,
    CriterionManager* const ptr_criterion_manager) {
    if (prog_name.empty()) {
        LogManager::LogWarning("Name of output file is not found, no flow field will be written in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        vtk_instance_.WriteVtuAll(prog_name, bool_binary_, k0OutputDataFormat_,
            *ptr_grid_manager, *ptr_criterion_manager);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
