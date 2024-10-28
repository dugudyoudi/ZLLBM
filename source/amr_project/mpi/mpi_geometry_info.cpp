//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
/**
 * @brief function to serializes a vector of geometry coordinates into a char buffer.
 * @param[in] vec_vertices vector store vertex information.
 * @param[out] ptr_buffer_size pointer to size of the buffer in bytes.
 * @return unique pointer to a char array to store the serialized data.
*/
std::unique_ptr<char[]> MpiManager::SerializeCoordiOrigin(
    const std::vector<std::unique_ptr<GeometryVertex>>& vec_vertices, int* const ptr_buffer_size) const {
    int num_points = 0;
    int real_size = sizeof(DefReal);
    if  (vec_vertices.size() * real_size* 3 > (std::numeric_limits<int>::max)()) {
        LogManager::LogError("size of the buffer is greater than the maximum of int"
            " in MpiManager::SerializeData(std::vector<GeometryCoordinate2D>) in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        num_points = static_cast<int>(vec_vertices.size());
    }
    int& buffer_size = *ptr_buffer_size;
    buffer_size = num_points * real_size * 3 + sizeof(int);
    // allocation buffer to store the serialized data
    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_points, sizeof(int));
    position += sizeof(int);
    for (auto& iter : vec_vertices) {
        std::memcpy(ptr_buffer + position, &iter->coordinate_.at(kXIndex), real_size);
        position += real_size;
        std::memcpy(ptr_buffer + position, &iter->coordinate_.at(kYIndex), real_size);
        position += real_size;
        std::memcpy(ptr_buffer + position, &iter->coordinate_.at(kZIndex), real_size);
        position += real_size;
    }
    return buffer;
}
/**
 * @brief function to deserializes data from a buffer into a vector of geometry coordinates.
 * @param[in] buffer unique pointer to a char array holding the serialized data.
 * @param[in] ptr_geo_info pointer to class storing geometry information.
 * @param[out] vec_points pointer to a vector store vertex information.
 */
void MpiManager::DeserializeCoordiOrigin(const std::unique_ptr<char[]>& buffer,
    GeometryInfoInterface* const ptr_geo_info,
    std::vector<std::unique_ptr<GeometryVertex>>* const ptr_vec_vertices) const {
    char* ptr_buffer = buffer.get();
    int real_size = sizeof(DefReal);
    int num_points;
    int position = 0;
    std::memcpy(&num_points, ptr_buffer, sizeof(int));
    position += sizeof(int);
    DefSizet size_pre = ptr_vec_vertices->size();
    ptr_vec_vertices->resize(size_pre + num_points);
    for (int i_node = 0; i_node < num_points; ++i_node) {
        ptr_vec_vertices->at(size_pre + i_node) = ptr_geo_info->GeoVertexCreator();
        std::memcpy(&(ptr_vec_vertices->at(size_pre + i_node)->coordinate_.at(kXIndex)),
            ptr_buffer + position, real_size);
        position += real_size;
        std::memcpy(&(ptr_vec_vertices->at(size_pre + i_node)->coordinate_.at(kYIndex)),
            ptr_buffer + position, real_size);
        position += real_size;
        std::memcpy(&(ptr_vec_vertices->at(size_pre + i_node)->coordinate_.at(kZIndex)),
            ptr_buffer + position, real_size);
        position += real_size;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI

