//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file io.h
* @author Zhengliang Liu
* @date  2022-8-14
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_IO_OUTPUT_FORMAT_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_IO_OUTPUT_FORMAT_H_
#include <string>
#include "../defs_libs.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace amrproject {
/**
* @class OutputDataFormatReal
* @brief class used to sepcify real format for output data.
* @note  change return data type of CastType for desired output format.
*/
class OutputDataFormatReal {
 public:
    std::string printf_format_{};
    std::string format_name_{};
    template<typename DataType>
    float CastType(const DataType& data) {
        return static_cast<float>(data);
    }
};
/**
* @class OutputDataFormatReal
* @brief class used to manage int format for output data.
*/
class OutputDataFormatInt {
 public:
    std::string printf_format_{};
    std::string format_name_{};
    template<typename DataType>
    std::int16_t CastType(const DataType& data) {
        return static_cast<int16_t>(data);
    }
};
/**
* @class OutputDataFormatSizet
* @brief class used to manage size_t format for output data.
*/
class OutputDataFormatSizet {
 public:
    std::string printf_format_{};
    std::string format_name_{};
    template<typename DataType>
    std::uint64_t CastType(const DataType& data) {
        return static_cast<uint64_t>(data);
    }
};
/**
* @class OutputDataFormatUint
* @brief class used to manage unsigned int format for output data.
*/
class OutputDataFormatUint {
 public:
    std::string printf_format_{};
    std::string format_name_{};
    template<typename DataType>
    std::uint16_t CastType(const DataType& data) {
        return static_cast<uint16_t>(data);
    }
};
/**
* @class OutputDataFormatReal
* @brief class used store output data.
*/
class OutputDataFormat {
public:
    OutputDataFormatReal output_real_;
    OutputDataFormatSizet output_sizet_;
    OutputDataFormatUint output_uint_;
    OutputDataFormatInt output_int_;
};
}  // end amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_IO_OUTPUT_FORMAT_H_
