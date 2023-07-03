//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file vtk_writer.h
* @author Zhengliang Liu
* @date  2022-8-04
* @brief  define the class to write vtk files.
*/

#ifndef ROOTPROJECT_SOURCE_IO_VTK_WRITER_H_
#define ROOTPROJECT_SOURCE_IO_VTK_WRITER_H_
#include <string>
#include <vector>
#include <unordered_map>
#include "../defs_libs.h"
#include "io/output_format.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
class Base64Utility {
 public:
    /**
    * @brief function to add data to end of a uint8 vector.
    * @param[in] data data need to convert to uint8 type.
    * @param[out] ptr_vec_uint8 a vector of uint8 type.
    */
    template <typename DataType>
    void AddToVecChar(const DataType& data,
        std::vector<uint8_t>* ptr_vec_uint8){
        const uint8_t* ptr_char = reinterpret_cast<const uint8_t*>(&data);
        ptr_vec_uint8->insert(ptr_vec_uint8->end(),
            ptr_char, ptr_char + sizeof(DataType));
    }
    void Encode(std::vector<uint8_t>* ptr_vec_uint8,
        std::vector<uint8_t>* ptr_vec_base64);

 private:
    static constexpr uint8_t kEncodeTable_[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
         "abcdefghijklmnopqrstuvwxyz0123456789+/";
};
/**
* @struct EVtkWriterGhostCellOption
* @brief enum class for methods to write ghots cells in vtk format
*/
enum class EVtkWriterGhostCellOption {
    kOutputEntirety = 0,  // grid at all given levels write in .vtu file
    kPartitionMultiBlock = 1, // grid write each block in .vtu file
    kPartitionEntirety = 2
};
/**
* @class VtkWriterManager
* @brief class used to manage vtk output.
* @date  2022-5-20
*/
class VtkWriterManager {
public:
    bool CheckIfLittleEndian() {
        short int test_ = 0x0001;
        return *((char*)&test_) ? true : false;
    };
    std::string str_vtk_byte_order_;
    std::string k0StrVtkAsciiOrBinary_;
    EVtkWriterGhostCellOption vtk_ghost_cell_option_ =
        EVtkWriterGhostCellOption::kOutputEntirety;

    void OptionInitial(const bool bool_binary);
    void WriteVtuAll(const std::string& program_name, const bool bool_binary,
        OutputDataFormat& output_data_format,
        GridManagerInterface& grid_manager,
        CriterionManager& criterion_manager);
    void WriteVtuGrid(const std::string& datafile_name,
        const bool bool_binary, const bool bool_overlap,
        const DefUint overlap_flag,
        const std::vector<DefSizet>& vec_level_in_one_vtu,
        OutputDataFormat& output_data_format,
        GridManagerInterface& grid_manager);
    void WriteVtuGeo(const std::string& datafile_name,
        const bool bool_binary,const bool bool_overlap,
        const DefUint overlap_flag, const DefUint dims,
        const std::vector<DefSizet>& vec_level_in_one_vtu,
        const std::array<DefReal,3>& space_offset,
        OutputDataFormat& output_data_format,
        CriterionManager& criterion_manager);

private:
    Base64Utility base64_instance_;
    const uint8_t kVtkNormalPoint_ = 0, kVtkDuplicatedPoint_ = 1,
        kVtkHiddenPoint_= 2;
  
    // definition from kitware
    // export const CellGhostTypes = {
    // DUPLICATECELL: 1, // the cell is present on multiple processors
    // HIGHCONNECTIVITYCELL : 2, // the cell has more neighbors than in a regular mesh
    // LOWCONNECTIVITYCELL : 4, // the cell has less neighbors than in a regular mesh
    // REFINEDCELL : 8, // other cells are present that refines it.
    // EXTERIORCELL : 16, // the cell is on the exterior of the data set
    // HIDDENCELL : 32, // the cell is needed to maintain connectivity, but the data values should be ignored.
    //  };

    //  export const PointGhostTypes = {
    //    DUPLICATEPOINT: 1, // the cell is present on multiple processors
    //    HIDDENPOINT : 2, // the point is needed to maintain connectivity, but the data values should be ignored.
    //  };

    void WritePvtu(FILE* const fp,
        const std::vector<std::string> vec_vtu_file_name,
        OutputDataFormat& output_data_format);

    void WriteGridPieces(FILE* const fp, const bool bool_binary,
        const bool bool_overlap, const DefUint overlap_flag,
        const GridInfoInterface& grid_info,
        OutputDataFormat& output_data_format,
        GridManagerInterface& grid_manager);
    void WirteGridCellOffsets(FILE* const fp, const bool bool_binary,
        const DefUint dims, const DefSizet num_cell,
        OutputDataFormat& output_data_format);
    void WirteGridCellTypes(FILE* const fp, const bool bool_binary,
        const DefUint dims, const DefSizet num_cell);
    void WirteGridNodeFlagStatus(FILE* const fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index,
        const DefMap<GridNode>& map_grid_node);
    void WirteGridNodeVtkVisualization(FILE* const fp, const bool bool_binary,
        const DefUint flag_overlap_visual,
        OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index,
        const DefMap<GridNode>& map_grid_node);

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    std::array<DefSizet, 2> CalculateNumOfGridCells(
        const bool bool_overlap, const DefUint overlap_flag,
        const SFBitsetAux2D& sfbitset_aux,
        const DefMap<GridNode>& map_grid_node,
        DefMap<DefSizet>* ptr_map_node_index) const;
    void WirteGridCoordinates(FILE* const fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const GridManager2D& grid_manager2d,
        const GridInfoInterface& grid_info,
        DefMap<DefSizet>* const ptr_map_node_index);
    void WirteGridCellConnectivity(FILE* const fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const SFBitsetAux2D& bitset_aux2d,
        const DefMap<DefSizet>& map_node_index);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
    void WirteGridCoordinates(FILE* const fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const GridManager3D& grid_manager3d,
        const GridInfoInterface& grid_info,
        DefMap<DefSizet>* const ptr_map_node_index);
    std::array<DefSizet, 2> CalculateNumOfGridCells(
        const bool bool_overlap, const DefUint overlap_flag,
        const SFBitsetAux3D& sfbitset_aux,
        const DefMap<GridNode>& map_grid_node,
        DefMap<DefSizet>* ptr_map_node_index) const;
    void WirteGridCellConnectivity(FILE* const fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const SFBitsetAux3D& bitset_aux3d,
        const DefMap<DefSizet>& map_node_index);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    void WriteGeometryPieces(FILE* const fp, const bool bool_binary,
        const DefUint dims, const std::array<DefReal, 3>& grid_offset,
        OutputDataFormat& output_data_format,
        GeometryInfoInterface& geo_info);
    DefSizet CalculateNumOfGeometryCells(
        const GeometryInfoInterface& geo_info) const;
    void WirteGeometryCellConnectivitykPolyLine(FILE* const fp,
        const bool bool_binary, const DefSizet num_points,
        OutputDataFormat& output_data_format);
    void WirteGeometryCellOffset(FILE* const fp, const bool bool_binary,
        const DefSizet num_cell, const EGeometryCellType geometry_cell_type,
        OutputDataFormat& output_data_format);
    void WirteGeometryCellType(FILE* const fp, const bool bool_binary,
        const DefSizet num_points, const EGeometryCellType geometry_cell_type);

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    DefSizet WirteGeometryCoordinates(FILE* const fp, const bool bool_binary,
        const DefUint dims, const std::array<DefReal, 3>& grid_offset,
        OutputDataFormat& output_data_format,
        const Geometry2DInterface& geo_info);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
    DefSizet WirteGeometryCoordinates(FILE* const fp, const bool bool_binary,
        const DefUint dims, const std::array<DefReal, 3>& grid_offset,
        OutputDataFormat& output_data_format,
        const Geometry3DInterface& geo_info);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
};
}  // end amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_IO_VTK_WRITER_MANAGER_H_
