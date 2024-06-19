//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
#include <limits>
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
#include "grid/grid_manager.h"
#ifdef DEBUG_UNIT_TEST
#include "../../googletest-main/googletest/include/gtest/gtest_prod.h"
#endif  // DEBUG_UNIT_TEST
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
    std::vector<DefSFCodeToUint> vec_sfcode_min_all_ranks_, vec_sfcode_max_all_ranks_;

    DefAmrIndexUint k0NumPartitionOuterLayers_ = 2;
    DefAmrIndexUint k0NumPartitionInnerLayers_ = k0NumPartitionOuterLayers_;
    std::vector<std::map<int, DefMap<DefAmrIndexUint>>> mpi_communication_inner_layers_;
    ///< inner layers (for sending) for mpi communication of all refinement levels
    /** nodes not on the partition interface and whose spacing filling code is
     between the minimum and the maximum spacing filling codes of current rank*/
    std::vector<DefMap<DefAmrIndexUint>> mpi_communication_outer_layers_;
    ///< outer layers (for receiving) for mpi communication of all refinement levels
    /** nodes whose spacing filling code exceeds 
     the minimum and the maximum spacing filling codes of current rank*/
    //  i      i       i  inner layer (rank 0)            o      o      o  outer layer (rank 1)
    //  -      -       -  inner layer: interface(rank 0)  o      o      o  outer layer (rank 1)
    //  o      o       o  outer layer (rank 0)            -      -      -  inner layer: interface (rank 1)
    //  o      o       o  outer layer (rank 0)            i      i      i  inner layer (rank 1)

    void SetUpMpi();
    void FinalizeMpi();

    void IniBroadcastBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds);
    inline int CheckNodeInWhichRank(const DefSFCodeToUint& sfcode_in) const {
        return static_cast<int>(std::lower_bound(vec_sfcode_max_all_ranks_.begin(),
            vec_sfcode_max_all_ranks_.end(), sfcode_in) - vec_sfcode_max_all_ranks_.begin());
    }

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
    std::unique_ptr<char[]> SerializeCoordiOrigin(
        const std::vector<GeometryCoordinate2D>& vec_points, int* const ptr_buffer_size) const;
    void DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate2D>* const vec_points) const;
    void IniSendNReceivePartitionedGeoCoordi(const std::array<DefReal, 2>& background_space,
        const SFBitsetAux2D& bitset_aux, const std::vector<DefSFBitset>& bitset_max,
        std::vector<GeometryCoordinate2D>* ptr_vec_coordinate);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
    std::unique_ptr<char[]> SerializeCoordiOrigin(
        const std::vector<GeometryCoordinate3D>& vec_points, int* const ptr_buffer_size) const;
    void DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate3D>* const vec_points) const;
    void IniSendNReceivePartitionedGeoCoordi(const std::array<DefReal, 3>& background_space,
        const SFBitsetAux3D& bitset_aux, const std::vector<DefSFBitset>& bitset_max,
        std::vector<GeometryCoordinate3D>* ptr_vec_coordinate);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    // functions to serialize and deserialize grid node information
 public:
    std::unique_ptr<char[]> SerializeNodeStoreUint(const DefMap<DefAmrUint>& map_nodes,
        int* const ptr_buffer_size) const;
    void DeserializeNodeStoreUint(const std::unique_ptr<char[]>& buffer,
        DefMap<DefAmrUint>* const map_nodes) const;
    std::unique_ptr<char[]> SerializeNodeSFBitset(const DefMap<DefAmrIndexUint>& map_nodes,
        int* const ptr_buffer_size) const;
    void DeserializeNodeSFBitset(const DefAmrIndexUint flag_node, const std::unique_ptr<char[]>& buffer,
        DefMap<DefAmrIndexUint>* const map_nodes) const;

 public:
    inline bool CodeLargerThanMax(const DefSFCodeToUint code_in, const DefSFCodeToUint code_max) const {
        return (code_in > code_max) ? true : false;
    }
    inline bool CodeSmallerThanMin(const DefSFCodeToUint code_in, const DefSFCodeToUint code_min) const {
        return (code_in < code_min) ? true : false;
    }
    // functions for inital partition
    void IniSendNReceiveTracking(const DefAmrIndexUint dims, const DefAmrIndexUint i_level,
        const std::vector<DefSFBitset>& bitset_max, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator,
        std::map<std::pair<ECriterionType, DefAmrIndexUint>, std::shared_ptr<TrackingGridInfoInterface>>*
        const ptr_map_tracking_info) const;
    void IniSendNReceiveCoarse2Fine0Interface(const DefAmrIndexUint dims, const DefAmrIndexUint i_level,
        const DefAmrIndexUint num_of_layers_coarse2fine, const DefAmrUint flag0,
        const DefMap<std::set<int>>& outmost_for_all_ranks,
        std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<InterfaceLayerInfo>>* const ptr_map_interface_info) const;
    void IniSendNReceivePartitionedGrid(const DefAmrIndexUint dims,
        const DefAmrIndexUint flag_size0, const DefAmrUint flag_coarse2fine,
        const std::vector<DefSFBitset>& bitset_min, const std::vector<DefSFBitset>& bitset_max,
        const std::vector<DefAmrIndexLUint>& indices_min, const std::vector<DefAmrIndexLUint>& indices_max,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset_rank0,
        std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_each,
        std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_ghost_each,
        std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const;
    void IniSendNReceiveGridInfoAtAllLevels(const DefAmrIndexUint flag_size0,
        const DefAmrUint flag_coarse2fine, const DefAmrIndexUint dims, const DefAmrIndexUint max_level,
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrIndexLUint>& indices_min, const std::vector<DefAmrIndexLUint>& indices_max,
        const std::vector<DefAmrIndexLUint>& vec_cost, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_creator,
        const std::vector<DefMap<DefAmrIndexUint>> ini_sfbitset_one_lower_level_rank0,
        std::array<DefSFBitset, 2>* const sfbitset_bound_current,
        std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_one_lower_level_current_rank,
        std::vector<DefMap<DefAmrIndexUint>>* const  ptr_sfbitset_ghost_one_lower_level_current_rank,
        std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info);

 private:
    void GridPartitionOnASingleRank(const DefMap<DefAmrIndexUint>& sfbitset_current_level);

// functions to serialize and deserialize information for node types other than grid node
 private:
    int CreateAndCommitCriterionIndexType(MPI_Datatype *ptr_mpi_pair_type) const;
    std::unique_ptr<char[]> IniSerializeTrackingNode(
        const std::set<DefSFCodeToUint>& set_nodes, int* const ptr_buffer_size) const;
    void IniDeserializeTrackingNode(const std::unique_ptr<char[]>& buffer,
        const TrackingNode& tracking_node_instance, DefMap<TrackingNode>* const ptr_map_tracking) const;
    std::unique_ptr<char[]> IniSerializeRefinementInterfaceNode(const DefMap<DefAmrUint>& interface_nodes,
        int* const ptr_buffer_size) const;
    void IniDeserializeRefinementInterfaceNode(const DefAmrUint flag0,
        const std::unique_ptr<char[]>& buffer, DefMap<DefAmrUint>* ptr_map_interface_layer) const;
    void IniSendNReiveOneLayerRefinementInterface(
        const DefAmrUint flag0, const DefMap<std::set<int>>& outmost_for_all_ranks,
        DefMap<DefAmrUint>* const ptr_map_interface_layer) const;

    // communicate between different partitioned blocks
 public:
    struct BufferSizeInfo {
        bool bool_exist_ = false;
        int num_chunks_ = 0;
        std::array<int, 2> array_buffer_size_ = {0, 0};
    };
    std::vector<bool> IdentifyRanksReceivingGridNode(const DefAmrIndexUint i_level) const;
    void SendNReceiveGridNodeBufferSize(const DefAmrIndexUint i_level,
        const int node_info_size, const std::vector<bool>& rank_id_sent,
        std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
        std::vector<BufferSizeInfo>* const ptr_receive_buffer_info) const;
    DefSizet CalculateBufferSizeForGridNode(const int rank_send,
        const GridInfoInterface& grid_inf) const;
    void SendGridNodeInformation(const int i_rank, const BufferSizeInfo& send_buffer_info,
        GridInfoInterface* ptr_grid_info, char* const ptr_buffer_send,
        std::vector<MPI_Request>* const ptr_reqs_send) const;
    std::unique_ptr<char[]> ReceiveGridNodeInformation(const int rank_receive, const int node_info_size,
        const BufferSizeInfo& receive_buffer_info, std::vector<MPI_Request>* const ptr_reqs_receive) const;
    void SendNReceiveGridNodes(
        std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
        std::vector<BufferSizeInfo>* const ptr_receive_buffer_info,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive,
        GridInfoInterface* ptr_grid_info) const;
    void WaitAndReadGridNodesFromBuffer(const std::vector<BufferSizeInfo>& send_buffer_info,
        const std::vector<BufferSizeInfo>& receive_buffer_info,
        const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        GridInfoInterface* ptr_grid_info) const;

 private:
    inline bool CheckBufferSizeNotExceedMax(DefSizet buffer_size) const {
        if (buffer_size > (std::numeric_limits<int>::max)()) {
            return false;
        } else {
            return true;
        }
    }

    // functions for the purpose of debug
 public:
    void CheckMpiNodesCorrespondence(const GridInfoInterface& grid_info) const;

 public:
#ifndef  DEBUG_DISABLE_2D_FUNCTION
    void IniTraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrIndexLUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
        const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void IniFindInterfaceForPartitionFromMinNMax(const DefSFCodeToUint& code_min,
        const DefSFCodeToUint& code_max, const std::array<DefAmrIndexLUint, 2>& code_domain_min,
        const std::array<DefAmrIndexLUint, 2>& code_domain_max, const SFBitsetAux2D& bitset_aux2d,
        DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const;
    void GetNLevelCorrespondingOnes2D(const DefAmrIndexUint i_level,
        const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_last_ones) const;
    int CheckNodeOnPartitionInterface2D(DefAmrIndexUint i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const DefSFBitset bitset_in, const SFBitsetAux2D& bitset_aux2d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background) const;
    void SearchForGhostLayerForMinNMax2D(const DefSFBitset bitset_in,
        const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_bound,
        bool (MpiManager::*ptr_func_compare)(const DefSFCodeToUint, const DefSFCodeToUint) const,
        const DefAmrIndexUint flag_ini, const SFBitsetAux2D& bitset_aux2d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        DefMap<DefAmrIndexUint>* const ptr_map_ghost_layer) const;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
    void IniTraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrIndexLUint>& vec_cost, const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
        const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void IniFindInterfaceForPartitionFromMinNMax(const DefSFCodeToUint& bitset_min,
        const DefSFCodeToUint& bitset_max, const std::array<DefAmrIndexLUint, 3>& code_domain_min,
        const std::array<DefAmrIndexLUint, 3>& code_domain_max, const SFBitsetAux3D& bitset_aux2d,
        DefMap<DefAmrIndexUint>* const ptr_partition_interface_background) const;
    void GetNLevelCorrespondingOnes3D(const DefAmrIndexUint i_level,
        const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_last_ones) const;
    int CheckNodeOnPartitionInterface3D(DefAmrIndexUint i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const DefSFBitset bitset_in, const SFBitsetAux3D& bitset_aux3d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background) const;
    void SearchForGhostLayerForMinNMax3D(const DefSFBitset bitset_in,
        const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_bound,
        bool (MpiManager::*ptr_func_compare)(const DefSFCodeToUint, const DefSFCodeToUint) const,
        const DefAmrIndexUint flag_ini, const SFBitsetAux3D& bitset_aux3d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        DefMap<DefAmrIndexUint>* const ptr_map_ghost_layer) const;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
};
}  //  end namespace amrproject
}  //  end namespace rootproject
#else  //  not define  ENABLE_MPI
namespace rootproject {
namespace amrproject {
/**
* @class MpiManager
* @brief class used to manage the mpi processes (empty).
*/
class MpiManager {
    // this is an empty class when mpi is not enabled.
    // used for argument passing in some functions.
};
}  //  end namespace amrproject
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_

