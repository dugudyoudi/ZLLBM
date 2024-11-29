//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_connection_to_grid.cpp
* @author Zhengliang Liu
* @brief functions to implement mpi communication for geometry information
* @date  2024-11-12
*/
#include "mpi/mpi_manager.h"
#include "criterion/geometry_info_interface.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to initialize vertices of the geometry on current rank.
* @param[in] code_min  minimum space filling code of current rank at background level.
* @param[in] code_max  maximum space filling code of current rank at background level.
* @param[in] sfbitset_aux   class manage space filling curves.
*/
void GeometryInfoInterface::InstantiateGeometryVertexInfo(
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max, const SFBitsetAuxInterface& sfbitset_aux) {
    std::vector<DefReal> grid_space_background = sfbitset_aux.GetBackgroundGridSpacing();
    std::vector<DefReal> coordi(k0GeoDim_, 0.);
    DefSFCodeToUint code_tmp;

    for (DefSizet i_vertex = 0; i_vertex < vec_vertices_.size(); ++i_vertex) {
        coordi[kXIndex] = vec_vertices_.at(i_vertex)->coordinate_[kXIndex];
        coordi[kYIndex] = vec_vertices_.at(i_vertex)->coordinate_[kYIndex];
        if (k0GeoDim_ == 3) {
            coordi[kZIndex] = vec_vertices_.at(i_vertex)->coordinate_[kZIndex];
        }
        code_tmp = sfbitset_aux.SFBitsetEncodingCoordi(grid_space_background, coordi).to_ullong();
        if (code_tmp >= code_min && code_tmp <= code_max) {
            map_vertices_info_.insert({ i_vertex, GeoInfoVertexCreator() });
            map_vertices_info_.at(i_vertex)->coordinate_.at(kXIndex) = coordi[kXIndex];
            map_vertices_info_.at(i_vertex)->coordinate_.at(kYIndex) = coordi[kYIndex];
            if (k0GeoDim_ == 3) {
                 map_vertices_info_.at(i_vertex)->coordinate_.at(kZIndex) = coordi[kZIndex];
            }
        }
    }
}
/**
 * @brief function to setup geometry for mpi communication.
 * @param[in] mpi_manager  class managing MPI communication.
 * @param[in] grid_info class storting grid node information on current rank.
 */
void GeometryInfoInterface::SetupGeometryInfo(DefReal time,
    const MpiManager& mpi_manager, const GridInfoInterface& grid_info) {
    InstantiateGeometryVertexInfo(mpi_manager.GetSFBitsetMinCurrentRank().to_ullong(),
        mpi_manager.GetSFBitsetMaxCurrentRank().to_ullong(), *grid_info.GetPtrSFBitsetAux());
    ptr_geo_shape_->UpdateShape(time);
}
}  // end namespace amrproject
}  // end namespace rootproject
