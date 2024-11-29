//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager.h
* @author Zhengliang Liu
* @date  2022-5-16
*/

#ifndef SOURCE_AMR_PROJECT_MPI_MPI_MANAGER_H_
#define SOURCE_AMR_PROJECT_MPI_MPI_MANAGER_H_
#include <bit>
#include <set>
#include <vector>
#include <functional>
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
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
#include "grid/grid_manager.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
* @class MpiManager
* @brief class used to manage the mpi processes.
*/
class MpiManager{
 protected:
    int num_of_ranks_ = 1;  ///< total number of mpi ranks
    int rank_id_;  ///< current rank
    DefSFBitset sfbitset_min_current_rank_, sfbitset_max_current_rank_;
    ///< space filling codes of background nodes on the interfaces of partitioned grid
    std::vector<DefSFCodeToUint> vec_sfcode_min_all_ranks_, vec_sfcode_max_all_ranks_;
    DefInt k0NumPartitionOuterLayers_ = 2;
    DefInt k0NumPartitionInnerLayers_ = k0NumPartitionOuterLayers_;

 public:
    /**
    * @class BufferSizeInfo
    * @brief class used to store information of buffer size.
    */  
    struct BufferSizeInfo {
        bool bool_exist_ = false;
        int num_chunks_ = 0;
        std::array<int, 2> array_buffer_size_ = {0, 0};
    };
    // set and get protected members
    int GetNumOfRanks() const { return num_of_ranks_; }
    int GetRankId() const { return rank_id_; }
    DefSFBitset GetSFBitsetMinCurrentRank() const { return sfbitset_min_current_rank_; }
    DefSFBitset GetSFBitsetMaxCurrentRank() const { return sfbitset_max_current_rank_; }
    std::vector<DefSFCodeToUint> GetSFCodeMinAllRanks() const { return vec_sfcode_min_all_ranks_; }
    std::vector<DefSFCodeToUint> GetSFCodeMaxAllRanks() const { return vec_sfcode_max_all_ranks_; }
    DefInt GetNumPartitionOuterLayers() const { return k0NumPartitionOuterLayers_; }
    DefInt GetNumPartitionInnerLayers() const { return k0NumPartitionInnerLayers_; }
    void SetSFBitsetMinCurrentRank(const DefSFBitset& sfbitset_min) { sfbitset_min_current_rank_ = sfbitset_min; }
    void SetSFBitsetMaxCurrentRank(const DefSFBitset& sfbitset_max) { sfbitset_max_current_rank_ = sfbitset_max; }
    void SetNumPartitionOuterLayers(const DefInt num_outer_layers) { k0NumPartitionOuterLayers_ = num_outer_layers; }
    void SetNumPartitionInnerLayers(const DefInt num_inner_layers) { k0NumPartitionInnerLayers_ = num_inner_layers; }
 
    std::vector<std::map<int, DefMap<DefInt>>> mpi_communication_inner_layers_;
    ///< inner layers (for sending) for mpi communication of all refinement levels
    /** nodes not on the partition interface and whose spacing filling code is
     between the minimum and the maximum spacing filling codes of current rank*/
    std::vector<DefMap<DefInt>> mpi_communication_outer_layers_;
    ///< outer layers (for receiving) for mpi communication of all refinement levels
    /** nodes whose spacing filling code exceeds 
     the minimum and the maximum spacing filling codes of current rank*/
    //  i      i       i  inner layer (rank 0)            o      o      o  outer layer (rank 1)
    //  -      -       -  inner layer: interface(rank 0)  o      o      o  outer layer (rank 1)
    //  o      o       o  outer layer (rank 0)            -      -      -  inner layer: interface (rank 1)
    //  o      o       o  outer layer (rank 0)            i      i      i  inner layer (rank 1)

    void SetUpMpi();
    void FinalizeMpi();

    void IniBroadcastSFBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds);
    /**
     * @brief function to find space filling code is on which rank.
     * @param[in] sfcode_in space filling code at background level.
     */
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
    std::unique_ptr<char[]> SerializeCoordiOrigin(
        const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices, int* const ptr_buffer_size) const;
    void DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer, GeometryInfoInterface* const ptr_geo_info,
        std::vector<std::unique_ptr<GeometryVertex>>* const ptr_vec_vertices) const;

    // functions to serialize and deserialize grid node information
 public:
    std::unique_ptr<char[]> SerializeNodeStoreInt(const DefMap<DefInt>& map_nodes,
        int* const ptr_buffer_size) const;
    void DeserializeNodeStoreInt(const std::unique_ptr<char[]>& buffer,
        DefMap<DefInt>* const map_nodes) const;
    std::unique_ptr<char[]> SerializeNodeSFBitset(const DefMap<DefInt>& map_nodes,
        int* const ptr_buffer_size) const;
    void DeserializeNodeSFBitset(const DefInt flag_node, const std::unique_ptr<char[]>& buffer,
        DefMap<DefInt>* const map_nodes) const;

 public:
    inline bool CodeLargerThanMax(const DefSFCodeToUint code_in, const DefSFCodeToUint code_max) const {
        return (code_in > code_max) ? true : false;
    }
    inline bool CodeSmallerThanMin(const DefSFCodeToUint code_in, const DefSFCodeToUint code_min) const {
        return (code_in < code_min) ? true : false;
    }
    // functions for inital partition
    void IniSendNReceiveTracking(const DefInt dims, const DefInt i_level,
        const std::vector<DefSFBitset>& bitset_max, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator,
        std::map<std::pair<ECriterionType, DefInt>, std::shared_ptr<TrackingGridInfoInterface>>*
        const ptr_map_tracking_info) const;
    void IniSendNReceiveCoarse2Fine0Interface(const DefInt dims, const DefInt i_level,
        const DefInt num_of_layers_coarse2fine, const DefInt flag0,
        const DefMap<std::set<int>>& outmost_for_all_ranks,
        std::map<std::pair<ECriterionType, DefInt>,
        std::shared_ptr<InterfaceLayerInfo>>* const ptr_map_interface_info) const;
    void IniSendNReceivePartitionedGrid(const DefInt dims,
        const DefInt flag_size0, const DefInt flag_coarse2fine,
        const std::vector<DefSFBitset>& bitset_min, const std::vector<DefSFBitset>& bitset_max,
        const std::vector<DefAmrLUint>& indices_min, const std::vector<DefAmrLUint>& indices_max,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefMap<DefInt>>& vec_sfbitset_rank0,
        std::vector<DefMap<DefInt>>* const ptr_sfbitset_each,
        std::vector<DefMap<DefInt>>* const ptr_sfbitset_c2f_each,
        std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const;
    void IniSendNReceiveGridInfoAtAllLevels(const DefInt flag_size0,
        const DefInt flag_coarse2fine, const DefInt dims, const DefInt max_level,
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefAmrLUint>& indices_min, const std::vector<DefAmrLUint>& indices_max,
        const std::vector<DefInt>& vec_cost, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_creator,
        const std::vector<DefMap<DefInt>> ini_sfbitset_one_lower_level_rank0,
        std::array<DefSFBitset, 2>* const sfbitset_bound_current,
        std::vector<DefMap<DefInt>>* const ptr_sfbitset_one_lower_level_current_rank,
        std::vector<DefMap<DefInt>>* const  ptr_sfbitset_ghost_one_lower_level_current_rank,
        std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info);

 private:
    void GridPartitionOnASingleRank(const DefMap<DefInt>& sfbitset_current_level);

// functions to serialize and deserialize information for node types other than grid node
 private:
    int CreateAndCommitCriterionIndexType(MPI_Datatype *ptr_mpi_pair_type) const;
    std::unique_ptr<char[]> IniSerializeTrackingNode(
        const std::set<DefSFCodeToUint>& set_nodes, int* const ptr_buffer_size) const;
    void IniDeserializeTrackingNode(const std::unique_ptr<char[]>& buffer,
        const TrackingNode& tracking_node_instance, DefMap<TrackingNode>* const ptr_map_tracking) const;
    std::unique_ptr<char[]> IniSerializeRefinementInterfaceNode(const DefMap<DefInt>& interface_nodes,
        int* const ptr_buffer_size) const;
    void IniDeserializeRefinementInterfaceNode(const DefInt flag0,
        const std::unique_ptr<char[]>& buffer, DefMap<DefInt>* ptr_map_interface_layer) const;
    void IniSendNReiveOneLayerRefinementInterface(
        const DefInt flag0, const DefMap<std::set<int>>& outmost_for_all_ranks,
        DefMap<DefInt>* const ptr_map_interface_layer) const;

    // communicate between different partitioned blocks
 public:
    void SendNReceiveSFbitsetForInterpolation(const DefInt i_level, const SFBitsetAuxInterface& sfbitset_aux,
        const DefMap<std::unique_ptr<GridNode>>& map_nodes_outer_layer,
        std::vector<int>* const ptr_vec_num_recv, std::map<int, DefMap<DefInt>>* const ptr_node_inner_layers) const;
    int SendNReceiveRequestNodesSFbitset(const std::vector<int>& num_node_request_current,
        const std::vector<int>& num_node_request_others, const std::vector<DefMap<DefInt>> &requested_nodes,
        std::map<int, DefMap<DefInt>>* const ptr_map_nodes_need_send) const;
    int SendGhostNodeForInterpolation(const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefInt>& vec_num_recv,
        const GridInfoInterface& grid_info_lower, GridInfoInterface* ptr_grid_info,
        std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
        std::vector<BufferSizeInfo>* const ptr_receive_buffer_info,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive) const;
    int WaitAndReadGhostNodeForInterpolation(const std::vector<BufferSizeInfo>& send_buffer_info,
        const std::vector<BufferSizeInfo>& receive_buffer_info,
        const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        GridInfoInterface* const ptr_grid_info) const;
    std::vector<bool> IdentifyRanksReceivingGridNode(const DefInt i_level) const;
    void SendNReceiveGridNodeBufferSize(const int node_info_size, const DefInt i_level,
        const std::map<int, DefMap<DefInt>>& map_send_nodes,
        std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
        std::vector<BufferSizeInfo>* const ptr_receive_buffer_info) const;
    DefSizet CalculateBufferSizeForGridNode(
        const DefSizet num_nodes, const GridInfoInterface& grid_info) const;
    void SendGridNodeInformation(const int i_rank, const BufferSizeInfo& send_buffer_info,
        const std::map<int, DefMap<DefInt>>& nodes_to_send,
        const std::function<void(const DefMap<DefInt>& , char* const)> func_copy_node_to_buffer,
        GridInfoInterface* ptr_grid_info, char* const ptr_buffer_send,
        std::vector<MPI_Request>* const ptr_reqs_send) const;
    std::unique_ptr<char[]> ReceiveGridNodeInformation(const int rank_receive, const int node_info_size,
        const BufferSizeInfo& receive_buffer_info, std::vector<MPI_Request>* const ptr_reqs_receive) const;
    void NonBlockingSendNReceiveGridNode(const int node_info_size,
        const std::map<int, DefMap<DefInt>>& map_send_nodes,
        const std::vector<BufferSizeInfo>& send_buffer_info,
        const std::vector<BufferSizeInfo>& receive_buffer_info,
        const std::function<void(const GridNode& node_ref, char* const)>& func_write_buffer,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive,
        GridInfoInterface* const ptr_grid_info) const;
    void SendAndReceiveGridNodesOnAllMpiLayers(std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
        std::vector<BufferSizeInfo>* const ptr_receive_buffer_info,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive,
        GridInfoInterface* ptr_grid_info) const;
    void WaitAndReadGridNodesFromBuffer(const std::vector<BufferSizeInfo>& send_buffer_info,
        const std::vector<BufferSizeInfo>& receive_buffer_info,
        const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
        const std::function<void(const char*,  GridNode* const)>& func_read_a_node_from_buffer,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
        std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
        GridInfoInterface* const ptr_grid_info) const;
    void MpiCommunicationForInterpolation(
        const SFBitsetAuxInterface& sfbitset_aux, const GridInfoInterface& grid_info_lower,
        GridInfoInterface* const ptr_grid_info) const;
    template <typename Node_T>
    std::unique_ptr<char[]> BlockingSendNReceiveGridNode(const int i_rank_send,
        const int i_rank_receive, const int node_info_size,
        const std::map<int, DefMap<DefInt>>& map_send_nodes,
        const std::vector<BufferSizeInfo>& send_buffer_info,
        const std::vector<BufferSizeInfo>& receive_buffer_info,
        const std::function<void(const Node_T& node_ref, char* const)>& func_write_buffer,
        const DefMap<std::unique_ptr<Node_T>>& map_nodes_info,
        std::vector<MPI_Request>* const ptr_vec_reqs_send,
        std::vector<std::unique_ptr<char[]>>* const ptr_vec_buffer_send)
        const requires std::is_base_of<GridNode, Node_T>::value;
    template <typename Node_T>
    void CopyNodeInfoToBuffer(const int node_info_size,
        const std::function<void(const Node_T& node_ref, char* const)>& func_copy_buffer,
        const DefMap<DefInt>& map_node_indices, const DefMap<std::unique_ptr<Node_T>>& map_nodes,
        char* const ptr_buffer) const requires std::is_base_of<GridNode, Node_T>::value;
    template <typename Node_T>
    void ReadNodeInfoFromBuffer(const int node_info_size,
        const std::function<void(const char*,  Node_T* const ptr_node)>& func_read_buffer,
        const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer,
        DefMap<std::unique_ptr<Node_T>>* const ptr_map_nodes)
        const requires std::is_base_of<GridNode, Node_T>::value;

 private:
    inline bool CheckBufferSizeNotExceedMax(DefSizet buffer_size) const {
        if (buffer_size > (std::numeric_limits<int>::max)()) {
            return false;
        } else {
            return true;
        }
    }

    // functions for the purpose of initialization
 public:
#ifndef  DEBUG_DISABLE_2D_FUNCTION
    void IniTraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefInt>& vec_cost, const std::vector<DefMap<DefInt>>& vec_sfbitset,
        const SFBitsetAux2D& bitset_aux2d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void IniFindInterfaceForPartitionFromMinNMax(const DefSFCodeToUint& code_min,
        const DefSFCodeToUint& code_max, const std::array<DefAmrLUint, 2>& code_domain_min,
        const std::array<DefAmrLUint, 2>& code_domain_max, const SFBitsetAux2D& bitset_aux2d,
        DefMap<DefInt>* const ptr_partition_interface_background) const;
    int CheckNodeOnPartitionInterface2D(DefInt i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const DefSFBitset bitset_in, const SFBitsetAux2D& bitset_aux2d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefInt>& partitioned_interface_background) const;
    void SearchForGhostLayerForMinNMax2D(const DefSFBitset bitset_in,
        const DefInt num_of_ghost_layers, const DefSFCodeToUint code_bound,
        bool (MpiManager::*ptr_func_compare)(const DefSFCodeToUint, const DefSFCodeToUint) const,
        const DefInt flag_ini, const SFBitsetAux2D& bitset_aux2d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        DefMap<DefInt>* const ptr_map_ghost_layer) const;
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
    void IniTraverseBackgroundForPartitionRank0(
        const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
        const std::vector<DefInt>& vec_cost, const std::vector<DefMap<DefInt>>& vec_sfbitset,
        const SFBitsetAux3D& bitset_aux3d, std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const;
    void IniFindInterfaceForPartitionFromMinNMax(const DefSFCodeToUint& bitset_min,
        const DefSFCodeToUint& bitset_max, const std::array<DefAmrLUint, 3>& code_domain_min,
        const std::array<DefAmrLUint, 3>& code_domain_max, const SFBitsetAux3D& bitset_aux2d,
        DefMap<DefInt>* const ptr_partition_interface_background) const;
    int CheckNodeOnPartitionInterface3D(DefInt i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const DefSFBitset bitset_in, const SFBitsetAux3D& bitset_aux3d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefInt>& partitioned_interface_background) const;
    void SearchForGhostLayerForMinNMax3D(const DefSFBitset bitset_in,
        const DefInt num_of_ghost_layers, const DefSFCodeToUint code_bound,
        bool (MpiManager::*ptr_func_compare)(const DefSFCodeToUint, const DefSFCodeToUint) const,
        const DefInt flag_ini, const SFBitsetAux3D& bitset_aux3d,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        DefMap<DefInt>* const ptr_map_ghost_layer) const;
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    // functions for the purpose of debug
 public:
    void DebugMpiForAllGrids(const GridManagerInterface& grid_manager) const;
    void CheckMpiNodesCorrespondence(const GridInfoInterface& grid_info) const;
    void CheckMpiPeriodicCorrespondence(const GridInfoInterface& grid_info) const;
};
/**
 * @brief function to copy node information to a buffer.
 * @param[in] node_info_size size of node information.
 * @param[in] func_copy_buffer function to copy node specified information to the buffer.
 * @param[in] map_node_indices container storing space filling codes of the nodes need to be copied.
 * @param[in] map_nodes_info container storing node information at a refinement level.
 * @param[out] ptr_buffer pointer to the buffer storing node information.
 */
template <typename Node_T>
void MpiManager::CopyNodeInfoToBuffer(const int node_info_size,
    const std::function<void(const Node_T& node_ref, char* const)>& func_copy_buffer,
    const DefMap<DefInt>& map_node_indices, const DefMap<std::unique_ptr<Node_T>>& map_nodes_info,
    char* const ptr_buffer) const requires std::is_base_of<GridNode, Node_T>::value {
    DefSizet position = 0;
    int key_size = sizeof(DefSFBitset);
    DefSizet buffer_size = (node_info_size + sizeof(DefSFBitset)) * map_node_indices.size();
    for (const auto& iter : map_node_indices) {
        if (map_nodes_info.find(iter.first) != map_nodes_info.end()) {
            std::memcpy(ptr_buffer + position, &(iter.first), key_size);
            position+=sizeof(DefSFBitset);
            func_copy_buffer(*map_nodes_info.at(iter.first), ptr_buffer + position);
            position+=node_info_size;
            if (position > buffer_size) {
                LogManager::LogError("Buffer to store node information overflows (buffer size is "
                    + std::to_string(buffer_size) + " ), please check"
                    " size of node info for mpi communication");
            }
        } else {
            LogManager::LogError("node with space filling code " + iter.first.to_string() + " does not exist");
        }
    }
}
/**
 * @brief function to read node information from a buffer.
 * @param[in] node_info_size size of node information.
 * @param[in] func_read_buffer function to read node specified information to the buffer.
 * @param[in] buffer_size total size of the buffer.
 * @param[in] buffer  pointer to the buffer storing node information.
 * @param[out] ptr_map_nodes pointer to container storing node information at a refinement level.
 */
template <typename Node_T>
void MpiManager::ReadNodeInfoFromBuffer(const int node_info_size,
    const std::function<void(const char*,  Node_T* const ptr_node)>& func_read_buffer,
    const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer,
    DefMap<std::unique_ptr<Node_T>>* const ptr_map_nodes)
    const requires std::is_base_of<GridNode, Node_T>::value {
    const DefSizet key_size = sizeof(DefSFBitset);
    const DefSizet num_nodes = buffer_size/(sizeof(DefSFBitset) + node_info_size);
    DefSizet position = 0;
    DefSFBitset key_code;
    const char* ptr_buffer = buffer.get();
    for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        if (ptr_map_nodes->find(key_code) != ptr_map_nodes->end()) {
            func_read_buffer(ptr_buffer + position, ptr_map_nodes->at(key_code).get());
        }   // may receive nodes that do not exist in current rank on c2f interface
        position += node_info_size;
        if (position > buffer_size) {
            LogManager::LogError("Buffer to store node information overflows, please check"
                " size of node info for mpi communication");
        }
    }
}
/**
 * @brief function to send and receive grid nodes via blocking mpi communication.
 * @param[in] node_info_size size of each node information will be sent.
 * @param[in] map_send_nodes  spacing filling codes of nodes will be sent.
 * @param[in] send_buffer_info buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] func_write_buffer function to write grid node infomation to buffer.
 * @param[in] map_nodes_info container storing node information at a refinement level.
 * @param[out] ptr_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_buffer_send pointer to buffer storing grid node information for sending.
 */
template <typename Node_T>
std::unique_ptr<char[]> MpiManager::BlockingSendNReceiveGridNode(const int i_rank_send,
    const int i_rank_receive, const int node_info_size,
    const std::map<int, DefMap<DefInt>>& map_send_nodes,
    const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::function<void(const Node_T&, char* const)>& func_write_buffer,
    const DefMap<std::unique_ptr<Node_T>>& map_nodes_info,
    std::vector<MPI_Request>* const ptr_vec_reqs_send,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_buffer_send) const
    requires std::is_base_of<GridNode, Node_T>::value {
    // send node information
    std::function<void(const DefMap<DefInt>& , char* const)> func_copy_node_to_buffer =
        [this, node_info_size, func_write_buffer, &map_nodes_info](
        const DefMap<DefInt>& map_send, char* const ptr_buffer) {
        CopyNodeInfoToBuffer<Node_T>(node_info_size, func_write_buffer, map_send, map_nodes_info, ptr_buffer);
    };

    int node_buffer_size = sizeof(DefSFBitset) + node_info_size;;
    if (send_buffer_info.at(i_rank_send).bool_exist_) {
        // copy grid node information to buffer
        ptr_vec_buffer_send->emplace_back(std::make_unique<char[]>(
            node_buffer_size*map_send_nodes.at(i_rank_send).size()));
        CopyNodeInfoToBuffer<Node_T>(node_info_size, func_write_buffer,
            map_send_nodes.at(i_rank_send), map_nodes_info, ptr_vec_buffer_send->back().get());

        // send grid node information chunk by chunk via non-blocking communication
        const int& buffer_size_send = send_buffer_info.at(i_rank_send).array_buffer_size_.at(0);
        const int& num_chunks = send_buffer_info.at(i_rank_send).num_chunks_;
        if (num_chunks > 1) {
            LogManager::LogWarning("buffer is not enough to send all grid node information at once on rank "
                + std::to_string(i_rank_send) + ", considering non-blocking communication.");
        }
        ptr_vec_reqs_send->resize(num_chunks);
        for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
            MPI_Isend(ptr_vec_buffer_send->back().get()+i_chunk*buffer_size_send,
                send_buffer_info.at(i_rank_send).array_buffer_size_.at(0),
                MPI_BYTE, i_rank_send, i_chunk, MPI_COMM_WORLD, &ptr_vec_reqs_send->at(i_chunk));
        }
        int i_chunk_last = num_chunks - 1;
        MPI_Isend(ptr_vec_buffer_send->back().get()+i_chunk_last*buffer_size_send,
            send_buffer_info.at(i_rank_send).array_buffer_size_.at(1), MPI_BYTE, i_rank_send, i_chunk_last,
            MPI_COMM_WORLD, &ptr_vec_reqs_send->at(i_chunk_last));
    }

    if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
        const int& num_chunks = receive_buffer_info.at(i_rank_receive).num_chunks_;
        const int& buffer_size_rest = receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(0);
        DefSizet buffer_size_total = (num_chunks - 1)*buffer_size_rest
            + receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1);
        std::unique_ptr<char[]> buffer_receive = std::make_unique<char[]>(buffer_size_total);
        int position = 0;
        for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
            MPI_Recv(buffer_receive.get()+i_chunk*buffer_size_rest, buffer_size_rest, MPI_BYTE, i_rank_receive,
                i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        int i_chunk_last = num_chunks - 1;
        MPI_Recv(buffer_receive.get()+ (num_chunks - 1)*buffer_size_rest,
            receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1), MPI_BYTE, i_rank_receive,
            i_chunk_last, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return buffer_receive;
    } else {
        return nullptr;
    }
}
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
 private:
    DefSFBitset sfbitset_min_current_rank_, sfbitset_max_current_rank_;

 public:
    DefSFBitset GetSFBitsetMinCurrentRank() const { return sfbitset_min_current_rank_; }
    DefSFBitset GetSFBitsetMaxCurrentRank() const { return sfbitset_max_current_rank_; }
    void SetSFBitsetMinCurrentRank(const DefSFBitset& sfbitset_min) { sfbitset_min_current_rank_ = sfbitset_min; }
    void SetSFBitsetMaxCurrentRank(const DefSFBitset& sfbitset_max) { sfbitset_max_current_rank_ = sfbitset_max; }
};
}  //  end namespace amrproject
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  SOURCE_AMR_PROJECT_MPI_MPI_MANAGER_H_

