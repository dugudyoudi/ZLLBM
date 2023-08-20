//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_geometry_info.cpp
* @author Zhengliang Liu
* @brief functions to communicate geometry information trough MPI
* @date  2022-8-5
*/
#include <limits>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
namespace rootproject {
namespace amrproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
 * @brief function to serializes a vector of GeometryCoordinate2D objects into a char buffer.
 * @param vec_points vector of GeometryCoordinate2D objects to be serialized.
 * @param buffer unique pointer to a char array to store the serialized data.
 * @return size of the buffer in bytes.
*/
int MpiManager::SerializeCoordiOrigin(const std::vector<GeometryCoordinate2D>& vec_points,
    std::unique_ptr<char[]>& buffer) const {
    int num_points;
    int real_size = sizeof(DefReal);
    if  (vec_points.size() * real_size* 2 > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater than the maximum of int"
            " in MpiManager::SerializeData(std::vector<GeometryCoordinate2D>)");
    } else {
        num_points = static_cast<int>(vec_points.size());
    }
    int buffer_size = num_points * real_size * 2 + sizeof(int);
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_points, sizeof(int));
    position += sizeof(int);
    for (auto& iter : vec_points) {
        std::memcpy(ptr_buffer + position, &iter.coordinate.at(kXIndex), real_size);
        position += real_size;
        std::memcpy(ptr_buffer + position, &iter.coordinate.at(kYIndex), real_size);
        position += real_size;
    }
    return buffer_size;
}
/**
 * @brief function to deserializes data from a buffer into a vector of 2D geometry coordinates.
 * @param[in] buffer unique pointer to a char array holding the serialized data.
 * @param[out] vec_points pointer to a vector store 2D geometry coordinates.
 */
void MpiManager::DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate2D>* const vec_points) const {
    char* ptr_buffer = buffer.get();
    int real_size = sizeof(DefReal);
    int num_points;
    int position = 0;
    std::memcpy(&num_points, ptr_buffer, sizeof(int));
    position += sizeof(int);
    DefSizet size_pre = vec_points->size();
    vec_points->resize(size_pre + num_points);
    for (int i_node = 0; i_node < num_points; ++i_node) {
        std::memcpy(&(vec_points->at(size_pre + i_node).coordinate.at(kXIndex)), ptr_buffer + position, real_size);
        position += real_size;
        std::memcpy(&(vec_points->at(size_pre + i_node).coordinate.at(kYIndex)), ptr_buffer + position, real_size);
        position += real_size;
    }
}
/**
 * @brief function to initializes the sending partitioned 2D geometry coordinates from rank 0 to others.
 * @param[in] background_space background grid space.
 * @param[in] bitset_aux class to handle functions for space filling curves.
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[out] ptr_vec_coordinate  pointer to a vector store 2D geometry coordinates.
 */
void MpiManager::IniSendNReceivePartitionedGeoCoordi(const std::array<DefReal, 2>& background_space,
    const SFBitsetAux2D& bitset_aux, const std::vector<DefSFBitset>& bitset_max,
    std::vector<GeometryCoordinate2D>* ptr_vec_coordinate) {
    DefSFBitset bitset_temp;
    std::array<DefAmrIndexLUint, 2> coordinate_index;
    const int rank_id = rank_id_, num_ranks = num_of_ranks_;
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefReal) / 2;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    DefSizet num_max = ull_max.size();
    GeometryCoordinate2D point_tmp;
    if (rank_id == 0) {
        std::vector<GeometryCoordinate2D> coordinate_origin_rank0;
        std::vector<std::vector<std::vector<GeometryCoordinate2D>>> vec_points_ranks(num_ranks);
        int index;
        std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
        for (const auto& i_point : *ptr_vec_coordinate) {
            coordinate_index =
            { static_cast<DefAmrIndexLUint>(i_point.coordinate.at(kXIndex) / background_space[kXIndex] + kEps),
              static_cast<DefAmrIndexLUint>(i_point.coordinate.at(kYIndex) / background_space[kYIndex] + kEps)};
            bitset_temp = bitset_aux.SFBitsetEncoding(coordinate_index);
            auto iter_index = std::lower_bound(ull_max.begin(),
            ull_max.end(), bitset_temp.to_ullong());
            index = static_cast<int>(iter_index - ull_max.begin());
#ifdef DEBUG_CHECK_GRID
            if (index == num_max) {
                LogManager::LogError("geometry point is out of computational "
                     "domain in GeometryInfo2DInterface::IniSendNReceivePartitionedGeoCoordiCoordi");
            }
#endif  // DEBUG_CHECK_GRID
            if (index > 0) {
                point_tmp.coordinate = {i_point.coordinate.at(kXIndex), i_point.coordinate.at(kYIndex)};
                if (i_counts.at(index) == 0) {
                    vec_points_ranks.at(index).push_back({point_tmp});
                    i_chunk_each_rank.at(index) += 1;
                    ++i_counts.at(index);
                } else {
                    vec_points_ranks.at(index).at(i_chunk_each_rank.at(index)).push_back(point_tmp);
                    ++i_counts.at(index);
                    if (i_counts.at(index) == max_buffer) {
                        // check if size of send buffer exceeds limits of int
                        i_counts.at(index) = 0;
                    }
                }
            } else {  // points on rank 0
                coordinate_origin_rank0.push_back(i_point);
            }
        }
        if (num_ranks > 1) {
            int buffer_size_send = 0;
            for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
                int num_chunks = static_cast<int>(vec_points_ranks.at(iter_rank).size());
                MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, 0, MPI_COMM_WORLD);
                std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                std::vector<MPI_Request> reqs_send(num_chunks);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    buffer_size_send = SerializeCoordiOrigin(
                           vec_points_ranks.at(iter_rank).at(i_chunk), vec_ptr_buffer.at(i_chunk));
                    MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                    i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
            }
            // delete coordinate_origin_ and initialize with partition one
            ptr_vec_coordinate->clear();
            *ptr_vec_coordinate = coordinate_origin_rank0;
        }
    } else {
        int num_chunks;
        MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
            int buffer_size_receive;
            MPI_Recv(&buffer_size_receive, 1, MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::unique_ptr<char[]> ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
            MPI_Recv(ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            DeserializeCoordiOrigin(ptr_buffer, ptr_vec_coordinate);
        }
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTION
/**
 * @brief function to serializes a vector of GeometryCoordinate2D objects into a char buffer.
 * @param vec_points vector of GeometryCoordinate2D objects to be serialized.
 * @param buffer unique pointer to a char array to store the serialized data.
 * @return size of the buffer in bytes.
*/
int MpiManager::SerializeCoordiOrigin(const std::vector<GeometryCoordinate3D>& vec_points,
    std::unique_ptr<char[]>& buffer) const {
    int num_points;
    int real_size = sizeof(DefReal);
    if  (vec_points.size() * real_size* 2 > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater than the maximum of int"
            " in MpiManager::SerializeData(std::vector<GeometryCoordinate3D>)");
    } else {
        num_points = static_cast<int>(vec_points.size());
    }
    int buffer_size = num_points * real_size * 3 + sizeof(int);
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_points, sizeof(int));
    position += sizeof(int);
    for (auto& iter : vec_points) {
        std::memcpy(ptr_buffer + position, &iter.coordinate.at(kXIndex), real_size);
        position += real_size;
        std::memcpy(ptr_buffer + position, &iter.coordinate.at(kYIndex), real_size);
        position += real_size;
        std::memcpy(ptr_buffer + position, &iter.coordinate.at(kZIndex), real_size);
        position += real_size;
    }
    return buffer_size;
}
/**
 * @brief function to deserializes data from a buffer into a vector of 2D geometry coordinates.
 * @param[in] buffer unique pointer to a char array holding the serialized data.
 * @param[out] vec_points pointer to a vector of GeometryCoordinate2D objects where the deserialized data will be stored.
 */
void MpiManager::DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
        std::vector<GeometryCoordinate3D>* const vec_points) const {
    char* ptr_buffer = buffer.get();
    int real_size = sizeof(DefReal);
    int num_points;
    int position = 0;
    std::memcpy(&num_points, ptr_buffer, sizeof(int));
    position += sizeof(int);
    DefSizet size_pre = vec_points->size();
    vec_points->resize(size_pre + num_points);
    for (int i_node = 0; i_node < num_points; ++i_node) {
        std::memcpy(&(vec_points->at(size_pre + i_node).coordinate.at(kXIndex)), ptr_buffer + position, real_size);
        position += real_size;
        std::memcpy(&(vec_points->at(size_pre + i_node).coordinate.at(kYIndex)), ptr_buffer + position, real_size);
        position += real_size;
        std::memcpy(&(vec_points->at(size_pre + i_node).coordinate.at(kZIndex)), ptr_buffer + position, real_size);
        position += real_size;
    }
}
/**
 * @brief function to initializes the sending partitioned 2D geometry coordinates from rank 0 to others.
 * @param[in] background_space background grid space.
 * @param[in] bitset_aux class to handle functions for space filling curves.
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[out] ptr_vec_coordinate  pointer to a vector store 2D geometry coordinates.
 */
void MpiManager::IniSendNReceivePartitionedGeoCoordi(const std::array<DefReal, 3>& background_space,
    const SFBitsetAux3D& sfbitset_aux, const std::vector<DefSFBitset>& bitset_max,
    std::vector<GeometryCoordinate3D>* ptr_vec_coordinate) {
    DefSFBitset bitset_temp;
    std::array<DefAmrIndexLUint, 3> coordinate_index;
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefReal) / 3;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    DefSizet num_max = ull_max.size();
    GeometryCoordinate3D point_tmp;
    if (rank_id == 0) {
        std::vector<GeometryCoordinate3D> coordinate_origin_rank0;
        std::vector<std::vector<std::vector<GeometryCoordinate3D>>> vec_points_ranks(num_ranks);
        int index;
        std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
        for (const auto& i_point : *ptr_vec_coordinate) {
            coordinate_index =
            { static_cast<DefAmrIndexLUint>(i_point.coordinate.at(kXIndex) / background_space[kXIndex] + kEps),
              static_cast<DefAmrIndexLUint>(i_point.coordinate.at(kYIndex) / background_space[kYIndex] + kEps),
              static_cast<DefAmrIndexLUint>(i_point.coordinate.at(kZIndex) / background_space[kZIndex] + kEps)};
            bitset_temp = sfbitset_aux.SFBitsetEncoding(coordinate_index);
            auto iter_index = std::lower_bound(ull_max.begin(),
            ull_max.end(), bitset_temp.to_ullong());
            index = static_cast<int>(iter_index - ull_max.begin());
#ifdef DEBUG_CHECK_GRID
            if (index == num_max) {
                LogManager::LogError("geometry point is out of computational "
                     "domain in MpiManager::IniSendNReceivePartitionedGeoCoordi");
            }
#endif  // DEBUG_CHECK_GRID
            if (index > 0) {
                point_tmp.coordinate = {i_point.coordinate.at(kXIndex),
                i_point.coordinate.at(kYIndex), i_point.coordinate.at(kZIndex)};
                if (i_counts.at(index) == 0) {
                    vec_points_ranks.at(index).push_back({point_tmp});
                    i_chunk_each_rank.at(index) += 1;
                    ++i_counts.at(index);
                } else {
                    vec_points_ranks.at(index).at(i_chunk_each_rank.at(index)).push_back(point_tmp);
                    ++i_counts.at(index);
                    if (i_counts.at(index) == max_buffer) {
                        // check if size of send buffer exceeds limits of int
                        i_counts.at(index) = 0;
                    }
                }
            } else {  // points on rank 0
                coordinate_origin_rank0.push_back(i_point);
            }
        }
        if (num_ranks > 1) {
            int buffer_size_send = 0;
            for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
                int num_chunks = static_cast<int>(vec_points_ranks.at(iter_rank).size());
                MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, 0, MPI_COMM_WORLD);
                std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                std::vector<MPI_Request> reqs_send(num_chunks);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    buffer_size_send = SerializeCoordiOrigin(
                           vec_points_ranks.at(iter_rank).at(i_chunk), vec_ptr_buffer.at(i_chunk));
                    MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                    i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
            }
            // delete coordinate_origin_ and initialize with partition one
            ptr_vec_coordinate->clear();
            *ptr_vec_coordinate = coordinate_origin_rank0;
        }
    } else {
        int num_chunks;
        MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
            int buffer_size_receive;
            MPI_Recv(&buffer_size_receive, 1, MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::unique_ptr<char[]> ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
            MPI_Recv(ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            DeserializeCoordiOrigin(ptr_buffer, ptr_vec_coordinate);
        }
    }
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI

