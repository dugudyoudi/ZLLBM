//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager.h
* @author Zhengliang Liu
* @date  2022-5-16
*/

#ifndef ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#define ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#include <bit>
#include <set>
#include <vector>
#include <map>
#include <utility>
#include <memory>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#ifdef __linux__
#include <netinet/in.h>
#include <endian.h>
#elif  _WIN32
#include <winsock2.h>
#endif
#include <mpi.h>
#include "criterion/geometry_coordi.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
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

    DefSFBitset sfbitset_min_current_rank_, sfbitset_max_current_rank_;
    ///< space filling codes of background nodes on the interfaces of partitioned grid

    void StartupMpi(int argc, char* argv[]);
    void FinalizeMpi();
    void SetMpiParameters();

    void IniBroadcastBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds);

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

    // criterion related functions
 public:
   
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    int SerializeCoordiOrigin(const std::vector<GeometryCoordinate2D>& vec_points,
        std::unique_ptr<char[]>& buffer) const;
    void DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate2D>* const vec_points) const;
    void IniSendNReceivePartitionedGeo(const std::array<DefReal, 2>& background_space,
        const SFBitsetAux2D& bitset_aux, const std::vector<DefSFBitset>& bitset_max,
        std::vector<GeometryCoordinate2D>* ptr_vec_coordinate);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
    int SerializeCoordiOrigin(const std::vector<GeometryCoordinate3D>& vec_points,
        std::unique_ptr<char[]>& buffer) const;
    void DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate3D>* const vec_points) const;
    void IniSendNReceivePartitionedGeo(const std::array<DefReal, 3>& background_space,
        const SFBitsetAux3D& bitset_aux, const std::vector<DefSFBitset>& bitset_max,
        std::vector<GeometryCoordinate3D>* ptr_vec_coordinate);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS


// grid related functions
 public:
    int SerializeNodeStoreUint(const DefMap<DefAmrUint>& map_nodes,
        std::unique_ptr<char[]>& buffer) const;
    void DeserializeNodeStoreUint(const std::unique_ptr<char[]>& buffer,
        DefMap<DefAmrUint>* const map_nodes) const;
    int SerializeNodeSFBitset(const DefMap<DefAmrIndexUint>& map_nodes,
        std::unique_ptr<char[]>& buffer) const;
    void DeserializeNodeSFBitset(const DefAmrIndexUint flag_node, const std::unique_ptr<char[]>& buffer,
        DefMap<DefAmrIndexUint>* const map_nodes) const;
    void IniSendNReceiveTracking(const DefAmrIndexUint dims, const DefAmrIndexUint i_level,
        const std::vector<DefSFBitset>& bitset_max, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator,
        std::map<std::pair<ECriterionType, DefAmrIndexUint>, std::shared_ptr<TrackingGridInfoInterface>>*
        const ptr_map_tracking_info) const;
    void IniSendNReceiveRefinementInterface(const DefAmrIndexUint dims, const DefAmrIndexUint i_level,
        const DefAmrIndexUint num_of_layers_coarse2fine, const std::vector<DefSFBitset>& bitset_max,
        const SFBitsetAuxInterface& sfbitset_aux, std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<InterfaceLayerInfo>>* const ptr_map_interface_info) const;
    void IniSendNReceivePartitionedGrid(const DefAmrIndexUint flag_size0,
        const std::vector<DefSFBitset>& bitset_max,
        const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_all_one_lower_level,
        const SFBitsetAuxInterface& sfbitset_aux,
        std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_each_one_lower_level) const;
    void SendNReceiveGridInfoAtGivenLevels(const DefAmrIndexUint flag_size0,
        const DefAmrIndexUint dims, const DefAmrIndexUint max_level,
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrIndexLUint>& vec_cost, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefMap<DefAmrIndexUint>>&  ini_sfbitset_one_lower_level_rank0,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_creator,
        std::array<DefSFBitset, 2>* const sfbitset_bound_current,
        std::vector<DefMap<DefAmrIndexUint>>* const  ptr_sfbitset_one_lower_level_current_rank,
        std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const;

 private:
    int CreateAndCommitCriterionIndexType(MPI_Datatype *ptr_mpi_pair_type) const;

    int IniSerializeTrackingNode(const std::set<DefSFCodeToUint>& set_nodes,
        std::unique_ptr<char[]>& buffer) const;
    void IniDeserializeTrackingNode(const std::unique_ptr<char[]>& buffer,
        const TrackingNode& tracking_node_instance, DefMap<TrackingNode>* const ptr_map_tracking) const;

    int IniSerializeRefinementInterfaceNode(const std::set<DefSFCodeToUint>& set_nodes,
        std::unique_ptr<char[]>& buffer) const;
    void IniDeserializeRefinementInterfaceNode(
        const std::unique_ptr<char[]>& buffer, DefMap<DefAmrUint>* ptr_map_interface_layer) const;
    void IniSendNReiveOneLayerRefinementInterface(const DefAmrIndexUint i_level,
        const std::vector<DefSFCodeToUint>& ull_max, const SFBitsetAuxInterface& sfbitset_aux,
        DefMap<DefAmrUint>* const ptr_map_interface_layer) const;

 public:
#ifndef  DEBUG_DISABLE_2D_FUNCTION
    void TraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
        const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void FindInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
        const DefSFBitset& bitset_max, const std::array<DefAmrIndexLUint, 2>& code_domain_min,
        const std::array<DefAmrIndexLUint, 2>& code_domain_max, const SFBitsetAux2D& bitset_aux2d,
        DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
    void TraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
        const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void FindInterfaceForPartitionFromMinNMax(const DefSFBitset& bitset_min,
        const DefSFBitset& bitset_max, const std::array<DefAmrIndexLUint, 3>& code_domain_min,
        const std::array<DefAmrIndexLUint, 3>& code_domain_max, const SFBitsetAux3D& bitset_aux2d,
        DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
};
}  //  end namespace amrproject
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_

