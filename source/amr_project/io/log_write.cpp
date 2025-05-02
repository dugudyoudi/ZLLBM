//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file log_write.cpp
* @author Zhengliang Liu
* @brief functions used to write log files.
* @date  2022-8-14
*/
#ifdef __cpp_lib_stacktrace
    #include <stacktrace>
    #define HAS_STACKTRACE 1
#else
    #define HAS_STACKTRACE 0
#endif
#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include <chrono>
#include <filesystem>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif  // ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
// static member
void LogManager::LogStartTime() {
    auto current_time = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());

    char time_char[26];
#ifdef _MSC_VER
    ctime_s(time_char, sizeof time_char, &current_time);
#else
    ctime_r(&current_time, time_char);
#endif

    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI

    FILE* fp = fopen((logfile_name_ + std::to_string(rank_id)).c_str(), "w");
    if (!fp) {
        printf("The file was not opened for writing log\n");
    } else {
        fprintf(fp, "project setup on: %s", time_char);
        if (rank_id == 0) {
            printf("project setup on: %s", time_char);
        }
#ifdef ENABLE_MPI
        int numprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        if (rank_id == 0) {
            printf("number of MPI ranks is: %d\n", numprocs);
        }
        fprintf(fp, "number of MPI ranks is: %d; current node is: %d \n", numprocs, rank_id);
#endif  // ENABLE_MPI
        fclose(fp);
    }
}
void LogManager::LogEndTime() {
    auto current_time = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());

    char time_char[26];

#ifdef _MSC_VER
    ctime_s(time_char, sizeof time_char, &current_time);
#else
    ctime_r(&current_time, time_char);
#endif

    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI

    FILE* fp = fopen((logfile_name_ + std::to_string(rank_id)).c_str(), "w");
    if (!fp) {
        printf("The file was not opened for writing log\n");
    } else {
        fprintf(fp, "project end on: %s", time_char);
        if (rank_id == 0) {
            printf("project end on: %s", time_char);
        }
        fclose(fp);
    }
}
/**
* @brief function to write normal information to the logfile.
* @param[in]  msg        information write to the logfile.
*/
void LogManager::LogInfo(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    FILE* fp = fopen((logfile_name_ + std::to_string(rank_id)).c_str(), "a");
    if (!fp) {
        printf("The log file was not opened\n");
    } else {
        fprintf(fp, "%s. \n", msg.c_str());
        fclose(fp);
    }
    printf("%s \n", msg.c_str());
}

/**
* @brief function to write warning information to the logfile.
* @param[in]  msg        information write to the logfile.
*/
void LogManager::LogWarning(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    FILE* fp = fopen((logfile_name_ + std::to_string(rank_id)).c_str(), "a");
    if (!fp) {
        printf("The log file was not opened\n");
    } else {
        fprintf(fp, "Warning: %s. \n", msg.c_str());
        fclose(fp);
    }
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
#endif
    printf("Warning of rank %d: %s \n", rank_id, msg.c_str());
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED
        | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif
}

/**
* @brief function to write error information to the logfile and stop the program.
* @param[in]  msg        information write to the logfile.
*/
void LogManager::LogError(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    FILE* fp = fopen((logfile_name_ + std::to_string(rank_id)).c_str(), "a");
#if HAS_STACKTRACE
    auto trace_msg = std::to_string(std::stacktrace::current());
    size_t first_newline_pos = trace_msg.find('\n');
    if (first_newline_pos != std::string::npos) {
        trace_msg = trace_msg.substr(first_newline_pos + 1);
    }
#else
    std::string trace_msg = "";
#endif

    if (!fp) {
        printf("The log file was not opened\n");
    } else {
        fprintf(fp, "Error: %s. \n%s", msg.c_str(), trace_msg.c_str());
        fclose(fp);
    }
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED);
#endif
    printf("Error of rank %d: %s. \n%s", rank_id, msg.c_str(), trace_msg.c_str());
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED
        | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif

#ifdef ENABLE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
}
}  // end namespace amrproject
}  // end namespace rootproject
