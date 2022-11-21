//  Copyright (c) 2022, Zhengliang Liu
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
#include "../defs_libs.h"
#include "io/output_format.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace amrproject {
class Base64Utility {
 public:
    template <typename DataType>
    void add_to_vec_char(const DataType& data,
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

    void OptionInitial(const bool bool_binary);
    void WriteVtu(const bool bool_binary,
        OutputDataFormat& output_data_format,
        std::shared_ptr<GridManagerInterface> ptr_grid_manager,
        std::shared_ptr<CriterionManager> ptr_criterion_manager);

private:
    Base64Utility base64_instance_;

    void WriteGrid(FILE* fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        std::shared_ptr<GridManagerInterface> ptr_grid_manager);
    std::array<DefLUint, 2> CalculateNumOfGridCells(
        const DefUint dims, const DefMap<GridNode>& map_grid_node,
        DefMap<DefLUint>* ptr_map_node_index);
    void WirteGridCoordinates(FILE* fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        std::shared_ptr<GridManagerInterface> ptr_grid_manager,
        const GridInfoInterface& grid_info);
    void WirteGridCellConnectivity(FILE* fp, const bool bool_binary,
        const DefUint dims, OutputDataFormat& output_data_format,
        const DefMap<GridNode>& map_grid_node,
        const DefMap<DefLUint>& map_node_index);
    void WirteGridCellOffsets(FILE* fp, const bool bool_binary,
        const DefUint dims, const DefLUint num_cell,
        OutputDataFormat& output_data_format);
    void WirteGridCellTypes(FILE* fp, const bool bool_binary,
        const DefUint dims, const DefLUint num_cell);


    void WriteGeometry(FILE* fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        std::shared_ptr<CriterionManager> ptr_criterion_manager);
    DefSizet CalculateNumOfGeometryCells(
        const GeometryInfoInterface& geo_info);
    void WirteGeometryCoordinates(FILE* fp, const bool bool_binary,
        const DefUint dims, OutputDataFormat& output_data_format,
        const GeometryInfoInterface& geo_info);
    void WirteGeometryCellConnectivity(FILE* fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const GeometryInfoInterface& geo_info);
    void WirteGeometryCellOffsets(FILE* fp, const bool bool_binary,
        OutputDataFormat& output_data_format,
        const GeometryInfoInterface& geo_info);
    void WirteGeometryCellTypes(FILE* fp, const bool bool_binary,
        const GeometryInfoInterface& geo_info);

};
}  // end amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_IO_VTK_WRITER_MANAGER_H_
