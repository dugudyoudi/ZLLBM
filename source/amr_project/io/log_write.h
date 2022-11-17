//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file log_write.h
* @author Zhengliang Liu
* @date  2022-8-14
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_IO_LOG_WRITE_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_IO_LOG_WRITE_H_
#include <string>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <cstdint>
namespace rootproject {
namespace amrproject {
namespace io {
void LogInfo(const std::string& msg);
void LogWarning(const std::string& msg);
void LogError(const std::string& msg);
void LogStartTime();
}  // end namespace io
}  // end amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_IO_LOG_WRITE_H_