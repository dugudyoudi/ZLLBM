//  Copyright (c) 2021 - 2023, Zhengliang Liu
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
#include <winsock2.h>
#include <Windows.h>
#endif
#include <cstdint>
namespace rootproject {
namespace amrproject {
class LogManager {
 public:
    inline static std::string logfile_name_{"log_of_node"};
    static void LogInfo(const std::string& msg);
    static void LogWarning(const std::string& msg);
    static void LogError(const std::string& msg);
    static void LogStartTime();
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_IO_LOG_WRITE_H_