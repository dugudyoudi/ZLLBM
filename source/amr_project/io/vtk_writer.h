//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file vtk_writer.h
* @author Zhengliang Liu
* @date  2022-8-04
* @brief  define the class to write vtk files.
*/

#ifndef SOURCE_AMR_PROJECT_IO_VTK_WRITER_H_
#define SOURCE_AMR_PROJECT_IO_VTK_WRITER_H_
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <functional>
#include "../defs_libs.h"
#include "io/output_format.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief enum class to convert data to uint8 type.
*/
class Base64Utility {
 public:
    /**
    * @brief function to add data to end of a uint8 vector.
    * @param[in] data data need to convert to uint8 type.
    * @param[out] ptr_vec_uint8 a vector of uint8 type.
    */
    template <typename DataType>
    void AddToVecChar(const DataType& data,
        std::vector<uint8_t>* ptr_vec_uint8) const {
        const uint8_t* ptr_char = reinterpret_cast<const uint8_t*>(&data);
        ptr_vec_uint8->insert(ptr_vec_uint8->end(),
            ptr_char, ptr_char + sizeof(DataType));
    }
    void Encode(std::vector<uint8_t>* ptr_vec_uint8,
        std::vector<uint8_t>* ptr_vec_base64) const;

 private:
    static constexpr uint8_t kEncodeTable_[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
         "abcdefghijklmnopqrstuvwxyz0123456789+/";
};
/**
* @struct EVtkWriterGhostCellOption
* @brief enum class for methods to write ghost cells in vtk format
*/
enum class EVtkWriterGhostCellOption {
    kOutputEntirety = 0,  // grid at all given levels write in .vtu file
    kPartitionMultiBlock = 1,  // grid write each block in .vtu file
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
        int16_t test_ = 0x0001;
        char* test_ptr = reinterpret_cast<char*>(&test_);
        return *test_ptr ? true : false;
    }
    std::string str_vtk_byte_order_;
    std::string k0StrVtkAsciiOrBinary_;
    EVtkWriterGhostCellOption vtk_ghost_cell_option_ = EVtkWriterGhostCellOption::kPartitionMultiBlock;
    Base64Utility base64_instance_;

    void OptionInitial(const bool bool_binary);
    void WriteVtuAll(const std::string& program_name, const bool bool_binary,
        const OutputDataFormat& output_data_format,
        const GridManagerInterface& grid_manager,
        const CriterionManager& criterion_manager);
    void WriteVtuGrid(const std::string& datafile_name,
        const bool bool_binary, const bool bool_overlap,
        const DefInt overlap_flag,
        const std::vector<DefInt>& vec_level_in_one_vtu,
        const OutputDataFormat& output_data_format,
        const GridManagerInterface& grid_manager);
    void WriteVtuGeo(const std::string& datafile_name,
        const bool bool_binary, const bool bool_overlap,
        const DefInt overlap_flag, const DefInt dims,
        const std::vector<DefInt>& vec_level_in_one_vtu,
        const std::array<DefReal, 3>& space_offset,
        const OutputDataFormat& output_data_format,
        const CriterionManager& criterion_manager);

 private:
    const uint8_t kVtkNormalPoint_ = 0, kVtkDuplicatedPoint_ = 1;

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

    void WritePvtu(FILE* const fp, const std::vector<std::string> vec_vtu_file_name,
        const OutputDataFormat& output_data_format,
        const std::vector<std::unique_ptr<OutputNodeVariableInfoInterface>>& output_variables);

    void WriteGridPieces(FILE* const fp, const bool bool_binary,
        const bool bool_overlap, const DefInt overlap_flag,
        const GridInfoInterface& grid_info,
        const OutputDataFormat& output_data_format,
        const GridManagerInterface& grid_manager);
    void WriteGridCellOffsets(FILE* const fp, const bool bool_binary,
        const DefInt dims, const DefSizet num_cell,
        const OutputDataFormat& output_data_format);
    void WriteGridCellTypes(FILE* const fp, const bool bool_binary,
        const DefInt dims, const DefSizet num_cell);
    void WriteGridNodeFlagStatus(FILE* const fp, const bool bool_binary,
        const OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index,
        const DefMap<std::unique_ptr<GridNode>>& map_grid_node);
    void WriteGridNodeVtkVisualization(FILE* const fp, const bool bool_binary,
        const DefInt flag_overlap_visual,
        const OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index,
        const DefMap<std::unique_ptr<GridNode>>& map_grid_node);

#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    std::array<DefSizet, 2> CalculateNumOfGridCells(
        const bool bool_overlap, const DefInt overlap_flag,
        const SFBitsetAux2D& sfbitset_aux,
        const DefMap<std::unique_ptr<GridNode>>& map_grid_node,
        DefMap<DefSizet>* ptr_map_node_index) const;
    void WriteGridCoordinates(FILE* const fp, const bool bool_binary,
        const OutputDataFormat& output_data_format, const GridManager2D& grid_manager2d,
        const GridInfoInterface& grid_info, DefMap<DefSizet>* const ptr_map_node_index);
    void WriteGridCellConnectivity(FILE* const fp, const bool bool_binary,
        const OutputDataFormat& output_data_format,
        const SFBitsetAux2D& bitset_aux2d, const DefMap<DefSizet>& map_node_index);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
    void WriteGridCoordinates(FILE* const fp, const bool bool_binary,
        const OutputDataFormat& output_data_format, const GridManager3D& grid_manager3d,
        const GridInfoInterface& grid_info, DefMap<DefSizet>* const ptr_map_node_index);
    std::array<DefSizet, 2> CalculateNumOfGridCells(const bool bool_overlap, const DefInt overlap_flag,
        const SFBitsetAux3D& sfbitset_aux, const DefMap<std::unique_ptr<GridNode>>& map_grid_node,
        DefMap<DefSizet>* ptr_map_node_index) const;
    void WriteGridCellConnectivity(FILE* const fp, const bool bool_binary,
        const OutputDataFormat& output_data_format,
        const SFBitsetAux3D& bitset_aux3d, const DefMap<DefSizet>& map_node_index);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    void WriteGeometryPieces(FILE* const fp, const bool bool_binary,
        const DefInt dims, const std::array<DefReal, 3>& grid_offset,
        const OutputDataFormat& output_data_format, const GeometryInfoInterface& geo_info);
    DefSizet CalculateNumOfGeometryCells(
        const GeometryInfoInterface& geo_info) const;
    void WriteGeometryCellConnectivityPolyLine(FILE* const fp, const bool bool_binary,
        const DefSizet num_points, const OutputDataFormat& output_data_format);
    void WriteGeometryCellOffset(FILE* const fp, const bool bool_binary,
        const DefSizet num_cell, const EGeometryCellType geometry_cell_type,
        const OutputDataFormat& output_data_format);
    void WriteGeometryCellType(FILE* const fp, const bool bool_binary,
        const DefSizet num_points, const EGeometryCellType geometry_cell_type);
    DefSizet WriteGeometryCoordinates(FILE* const fp, const bool bool_binary,
        const DefInt dims, const std::array<DefReal, 3>& grid_offset,
        const OutputDataFormat& output_data_format, const GeometryInfoInterface& geo_info);

};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_IO_VTK_WRITER_H_
