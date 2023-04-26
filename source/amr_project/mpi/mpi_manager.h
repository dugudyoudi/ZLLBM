//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file mpi.h
* @author Zhengliang Liu
* @date  2022-5-16
*/

#ifndef ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#define ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#include <bit>
#include <vector>
#include <memory>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#ifdef __linux__
#include <netinet/in.h>
#include <endian.h>
#endif
#include <mpi.h>
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @class MpiManager
* @brief class used to manage the mpi processes.
*/
class MpiManager{
 public:
    int num_of_ranks_ = 1;  ///< total number of mpi ranks
    int rank_id_;  ///< current rank
    void StartupMpi(int argc, char* argv[]);
    void SetMpiParameters();

    int SerializeData(const DefMap<DefUint>& map_nodes,
        std::unique_ptr<char[]>& buffer) const;
    void DeserializeData(const std::unique_ptr<char[]>& buffer,
        DefMap<DefUint>* const map_nodes) const;

    void IniPartiteGridBySpaceFillingCurves(
        const std::vector<DefMap<DefUint>>& sfbitset_one_lower_level,
        GridManagerInterface const& grid_manager,
        std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max);
    void IniBroadcastBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds);
    void IniSendNReceivePartitionedGrid(
        const std::vector<DefSFBitset>& ptr_bitset_max,
        const std::vector<DefMap<DefUint>>& sfbitset_all_one_lower_level,
        GridManagerInterface const& grid_manager,
        std::vector<DefMap<DefUint>>* const ptr_sfbitset_each_one_lower_level) const;

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    void IniSendNReceivePartitionedGeoCoordi(
        const GridManager2D& grid_manager2d,
        const std::vector<DefSFBitset>& bitset_max,
        Geometry2DInterface* const ptr_geo2d) const;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
    void IniSendNReceivePartitionedGeoCoordi(
        const GridManager3D& grid_manager3d,
        const std::vector<DefSFBitset>& bitset_max,
        Geometry3DInterface* const ptr_geo3d) const;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    // functions to convert host unsigned integer to network one
    inline uint8_t HtoNUint(uint8_t val_host) const {
        return val_host;
    }
    inline uint16_t HtoNUint(uint16_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
            return htons(val_host);
        } else {
            return val_host;
        }
    }
    inline uint32_t HtoNUint(uint32_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
            return htonl(val_host);
        } else {
            return val_host;
        }
    }
    inline uint64_t HtoNUint(uint64_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
        return htonll(val_host);
#elif __linux__
        return htobe64(val_host);
#endif
        } else {
            return val_host;
        }
    }
    // functions to convert network unsigned integer to host one
    inline uint8_t NtoHUint(uint8_t val_network) const {
        return val_network;
    }
    inline uint16_t NtoHUint(uint16_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
            return ntohs(val_network);
        } else {
            return val_network;
        }
    }
    inline uint32_t NtoHUint(uint32_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
            return ntohl(val_network);
        } else {
            return val_network;
        }
    }
    inline uint64_t NtoHUint(uint64_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
        return ntohll(val_network);
#elif __linux__
        return be64toh(val_network);
#endif
        } else {
            return val_network;
        }
    }
    // functions to convert network integer to host one
    inline uint8_t HtoNInt(int8_t val_host) const {
        return *reinterpret_cast<uint8_t*>(&val_host);
    }
    inline uint16_t HtoNInt(int16_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
            return htons(*reinterpret_cast<uint16_t*>(&val_host));
        } else {
            return *reinterpret_cast<uint16_t*>(&val_host);
        }
    }
    inline uint32_t HtoNInt(int32_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
            return htonl(*reinterpret_cast<uint32_t*>(&val_host));
        } else {
            return *reinterpret_cast<uint32_t*>(&val_host);
        }
    }
    inline uint64_t HtoNInt(int64_t val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
        return htonll(*reinterpret_cast<uint64_t*>(&val_host));
#elif __linux__
        return htobe64(*reinterpret_cast<uint64_t*>(&val_host));
#endif
        } else {
            return *reinterpret_cast<uint64_t*>(&val_host);
        }
    }
    // functions to convert network integer to host one
    inline int8_t NtoHInt(uint8_t val_network) const {
        return *reinterpret_cast<int8_t*>(&val_network);
    }
    inline int16_t NtoHInt(uint16_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
            uint16_t val_tmp = ntohs(val_network);
            return *reinterpret_cast<int16_t*>(&val_tmp);
        } else {
            return *reinterpret_cast<int16_t*>(&val_network);
        }
    }
    inline int32_t NtoHInt(uint32_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
            uint32_t val_tmp = ntohl(val_network);
            return *reinterpret_cast<int32_t*>(&val_tmp);
        } else {
            return *reinterpret_cast<int32_t*>(&val_network);
        }
    }
    inline int64_t NtoHInt(uint64_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
            uint64_t val_tmp = ntohll(val_network);
            return *reinterpret_cast<int64_t*>(&val_tmp);
#elif __linux__
            uint64_t val_tmp = be64toh(val_network);
            return *reinterpret_cast<int64_t*>(&val_tmp);
#endif
        } else {
            return *reinterpret_cast<int64_t*>(&val_network);
        }
    }
    // functions to convert network real to host one
    inline uint32_t HtoNReal(float val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
            return htonl(*reinterpret_cast<uint32_t*>(&val_host));
        } else {
            return *reinterpret_cast<uint32_t*>(&val_host);
        }
    }
    inline uint64_t HtoNReal(double val_host) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
        return htonll(*reinterpret_cast<uint64_t*>(&val_host));
#elif __linux__
        return htobe64(*reinterpret_cast<uint64_t*>(&val_host));
#endif
        } else {
            return *reinterpret_cast<uint64_t*>(&val_host);
        }
    }
    // functions to convert host real to network one
    inline float NtoHReal(uint32_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
            uint32_t val_tmp = ntohl(val_network);
            return *reinterpret_cast<float*>(&val_tmp);
        } else {
            return *reinterpret_cast<float*>(&val_network);
        }
    }
    inline double NtoHReal(uint64_t val_network) const {
        if constexpr (std::endian::native == std::endian::little) {
#ifdef WIN32
            uint64_t val_tmp = ntohll(val_network);
            return *reinterpret_cast<double*>(&val_tmp);
#elif __linux__
            uint64_t val_tmp = be64toh(val_network);
            return *reinterpret_cast<double*>(&val_tmp);
#endif
        } else {
            return *reinterpret_cast<double*>(&val_network);
        }
    }
};
}  //  end namespace amrproject
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_

