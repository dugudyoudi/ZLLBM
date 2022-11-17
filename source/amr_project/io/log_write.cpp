//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file log_write.cpp
* @author Zhengliang Liu
* @date  2022-8-14
*/
#include <chrono>
#include "io/log_write.h"
#include "io/io_manager.h"
namespace rootproject {
namespace amrproject {
namespace io {
// static member
std::string IoManager::logfile_name;

void LogStartTime() {
    auto current_time = std::chrono::system_clock::to_time_t
    (std::chrono::system_clock::now());

    char time_char[26];
    ctime_s(time_char, sizeof time_char, &current_time);

    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::rank_id_;
#endif  // ENABLE_MPI

    FILE* fp = nullptr;
    errno_t err = fopen_s(&fp,
        (IoManager::logfile_name + std::to_string(rank_id)).c_str(), "w");
    if (!fp) {
        printf_s("The file was not opened for writing log\n");
    } else {
        fprintf_s(fp, "project setup on: %s", time_char);
        fclose(fp);
    }
    if (rank_id == 0) printf("project setup on: %s", time_char);
}
/**
* @brief function to write normal information to the logfile.
* @param[in]  msg        information write to the logfile.
*/
void LogInfo(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::rank_id_;
#endif  // ENABLE_MPI
    FILE* fp = nullptr;
    errno_t err = fopen_s(&fp,
        (IoManager::logfile_name + std::to_string(0)).c_str(), "a");
    if (!fp) {
        printf_s("The log file was not opened\n");
    } else {
        fprintf(fp, "%s \n", msg.c_str());
        fclose(fp);
    }
    printf("%s \n", msg.c_str());
}
/**
* @brief function to write normal information to the logfile.
* @param[in]  mpi_id        id of the mpi node.
* @param[in]  msg        information write to the logfile.
*/
void LogInfo(const int mpi_id, const std::string& msg) {
    FILE* fp = nullptr;
    errno_t err = fopen_s(&fp,
        (IoManager::logfile_name + std::to_string(mpi_id)).c_str(), "a");
    if (!fp) {
        printf_s("The log file was not opened\n");
    } else {
        fprintf(fp, "Infor: %s \n", msg.c_str());
        fclose(fp);
    }
    printf("Infor of node %d: %s \n", mpi_id, msg.c_str());
}

/**
* @brief function to write warning information to the logfile.
* @param[in]  mpi_id        id of the mpi node.
* @param[in]  msg        information write to the logfile.
*/
void LogWarning(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::rank_id_;
#endif  // ENABLE_MPI
    FILE* fp = nullptr;
    errno_t err = fopen_s(&fp,
        (IoManager::logfile_name + std::to_string(rank_id)).c_str(), "a");
    if (!fp) {
        printf_s("The log file was not opened\n");
    } else {
        fprintf(fp, "Warning: %s \n", msg.c_str());
        fclose(fp);
    }
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
#endif
    printf("Warning of node %d: %s \n", rank_id, msg.c_str());
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED
        | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif
}

/**
* @brief function to write error information to the logfile.
* @param[in]  mpi_id        id of the mpi node.
* @param[in]  msg        information write to the logfile.
*/
void LogError(const std::string& msg) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::rank_id_;
#endif  // ENABLE_MPI
    FILE* fp = nullptr;
    errno_t err = fopen_s(&fp,
        (IoManager::logfile_name + std::to_string(rank_id)).c_str(), "a");
    if (!fp) {
        printf_s("The log file was not opened\n");
    } else {
        fprintf(fp, "Error: %s \n",  msg.c_str());
        fclose(fp);
    }
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED);
#endif
    printf("Error of node %d: %s \n", rank_id, msg.c_str());
#ifdef _WIN32
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
        FOREGROUND_INTENSITY | FOREGROUND_RED
        | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif
    exit(0);
}
}  // end namespace io
}  // end amrproject
}  // end namespace rootproject
