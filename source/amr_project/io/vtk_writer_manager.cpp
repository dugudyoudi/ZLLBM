//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file vtk_writer_manager.cpp
* @author Zhengliang Liu
* @brief functions used to write vtk data.
* @date  2022-8-04
*/
#include <typeinfo>
#include <filesystem>
#include "auxiliary_inline_func.h"
#include "io/log_write.h"
#include "io/io_manager.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief  funciton to convert data to base64
         (split three bytes into four six-bit values)
* @param[in]  ptr_vec_uint8 input data
* @param[in]  ptr_vec_base64 converted base64 data
*/
void Base64Utility::Encode(std::vector<uint8_t>* ptr_vec_uint8,
    std::vector<uint8_t>* ptr_vec_base64) {
    // number of bytes of the output data
    const uint32_t size_vec = static_cast<uint32_t>(ptr_vec_uint8->size());
    ptr_vec_uint8->insert(ptr_vec_uint8->begin(),
        reinterpret_cast<const uint8_t*> (&size_vec),
        reinterpret_cast<const uint8_t*> (&size_vec) + sizeof(uint32_t));
    DefSizet uint8_length = ptr_vec_uint8->size();
    DefSizet base64_length = uint8_length / 3 * 4;
    if (uint8_length % 3 > 0) {
        base64_length += 4;
    }
    ptr_vec_base64->resize(base64_length);

    DefSizet i, index = 0;
    for (i = 0; i + 3 <= uint8_length; i += 3) {
        ptr_vec_base64->at(index++) =
            kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
        ptr_vec_base64->at(index++) =
            kEncodeTable_[((ptr_vec_uint8->at(i) << 4) & 0x30)
            | ((ptr_vec_uint8->at(i + 1) >> 4) & 0x0F)];
        ptr_vec_base64->at(index++) =
            kEncodeTable_[((ptr_vec_uint8->at(i + 1) << 2) & 0x3C)
            | (ptr_vec_uint8->at(i + 2) >> 6) & 0x03];
        ptr_vec_base64->at(index++) =
            kEncodeTable_[ptr_vec_uint8->at(i + 2) & 0x3F];
    }

    if (i < uint8_length) {
        DefSizet bytes_remain = uint8_length - i;
        if (bytes_remain == 1) {  // single byte
            ptr_vec_base64->at(index++) =
                kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
            ptr_vec_base64->at(index++) =
                kEncodeTable_[(ptr_vec_uint8->at(i) << 4) & 0x30];
            ptr_vec_base64->at(index++) = '=';
            ptr_vec_base64->at(index++) = '=';
        } else  {  // pair bytes
            ptr_vec_base64->at(index++) =
                kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
            ptr_vec_base64->at(index++) =
                kEncodeTable_[((ptr_vec_uint8->at(i) << 4) & 0x30)
                | ((ptr_vec_uint8->at(i + 1) >> 4)) & 0x0F];
            ptr_vec_base64->at(index++) =
                kEncodeTable_[(ptr_vec_uint8->at(i + 1) << 2) & 0x3C];
            ptr_vec_base64->at(index++) = '=';
        }
    }
}
/**
* @brief  funciton to initialize vtk related options
* @param[in]  bool_binary write data in bindary or ascii format
*/
void VtkWriterManager::OptionInitial(const bool bool_binary) {
    // machine dependent option
    if (CheckIfLittleEndian()) {
        str_vtk_byte_order_ = "LittleEndian";
    } else {
        str_vtk_byte_order_ = "BigEndian";
    }
    if (bool_binary) {
        k0StrVtkAsciiOrBinary_ = "binary";
    } else {
        k0StrVtkAsciiOrBinary_ = "ascii";
    }
}
/**
* @brief   function to write all grid levels
* @param[in]  folder_name name of the folder
* @param[in]  bool_binary write data in binary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  grid_manager class to manage grid information
* @param[in]  criterion_manager class to manage criterion information
*/
void VtkWriterManager::WriteVtuAll(const std::string& folder_name,
    const bool bool_binary, OutputDataFormat& output_data_format,
    GridManagerInterface& grid_manager,
    CriterionManager& criterion_manager) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    std::filesystem::create_directories(folder_name);
    std::string proj_and_rank = folder_name + "/rank_"
        + std::to_string(rank_id) + "_";

    std::array<DefReal, 3> grid_offset;
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
    if (grid_manager.k0GridDims_ == 2) {
        const GridManager2D& grid_manager2d =
            dynamic_cast<GridManager2D&> (grid_manager);
        grid_offset.at(kXIndex) =
            grid_manager2d.k0RealOffset_.at(kXIndex);
        grid_offset.at(kYIndex) =
            grid_manager2d.k0RealOffset_.at(kYIndex);
        grid_offset.at(kZIndex) = 0.;
    }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
    if (grid_manager.k0GridDims_ == 3) {
        const GridManager3D& grid_manager3d =
            dynamic_cast<GridManager3D&> (grid_manager);
        grid_offset.at(kXIndex) =
            grid_manager3d.k0RealOffset_.at(kXIndex);
        grid_offset.at(kYIndex) =
            grid_manager3d.k0RealOffset_.at(kYIndex);
        grid_offset.at(kZIndex) =
            grid_manager3d.k0RealOffset_.at(kZIndex);
    }
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    bool bool_overlap;
    DefAmrUint overlap_flag = 0;
    switch (vtk_ghost_cell_option_) {
    case EVtkWriterGhostCellOption::kOutputEntirety : {
        bool_overlap = false;
        overlap_flag = grid_manager.kNodeStatusFine2CoarseM1_
            | grid_manager.kNodeStatusFine2Coarse0_;
        std::vector<DefAmrIndexUint> output_levels;
        for (DefAmrIndexUint i = 0; i < grid_manager.k0MaxLevel_ + 1; ++i) {
            output_levels.push_back(i);
        }
        std::string grid_file_name = proj_and_rank + "gridmerged_"
            + std::to_string(1);
        WriteVtuGrid(grid_file_name, bool_binary, bool_overlap, overlap_flag,
            output_levels, output_data_format, grid_manager);
        break;
    }
    case EVtkWriterGhostCellOption::kPartitionMultiBlock: {
        std::vector<DefAmrIndexUint> output_levels;
        bool_overlap = true;
        overlap_flag = grid_manager.kNodeStatusCoarse2FineM1_
            | grid_manager.kNodeStatusCoarse2FineM1_;
        std::vector<std::string> vec_vtu_file_name;
        for (DefAmrIndexUint i = 0; i < grid_manager.k0MaxLevel_ + 1; ++i) {
            std::string grid_file_name = proj_and_rank + "grid_level_"
                + std::to_string(i);
            WriteVtuGrid(grid_file_name, bool_binary, bool_overlap, overlap_flag,
                {i}, output_data_format, grid_manager);
            vec_vtu_file_name.push_back(grid_file_name);
        }

        // write pv
        FILE* fp = nullptr;
        errno_t err = fopen_s(&fp, (folder_name + ".pvtu").c_str(), "w");
        if (!fp) {
            LogError("File was not opened for writing"
                " pvtu file in WriteVtuAll.");
        } else {
            WritePvtu(fp, vec_vtu_file_name, output_data_format);
        }
        break;
    }
    default:
        break;
    }

    bool_overlap = false;
    std::vector<DefAmrIndexUint> output_geos;
    for (DefAmrIndexUint i = 0; i < criterion_manager.vec_ptr_geometries_.size(); ++i) {
        output_geos.push_back(i);
    }
    std::string geo_file_name = proj_and_rank + "mergedgeo_"
        + std::to_string(1);
    WriteVtuGeo(geo_file_name, bool_binary, bool_overlap, overlap_flag,
        grid_manager.k0GridDims_, output_geos, grid_offset,
        output_data_format, criterion_manager);
}
/**
* @brief   function to write xml formatted unstructed mesh (.vtu)
* @param[in]  bool_binary write data in bindary or ascii format.
* @param[in]  bool_overlap write overlapping region or not.
* @param[in]  overlap_flag flag indicates nodes in the overlapping region.
* @param[in]  dims dimesion of mesh.
* @param[in]  vec_level_in_one_vtu idices of geometries write in this vtu file.
* @param[in]  grid_offset offset coordinates of the mesh.
* @param[in]  output_data_format output data (real or integer) format.
* @param[in]  grid_manager class to manage grid information.
* @param[in]  criterion_manager class to manage criterion information.
*/
void VtkWriterManager::WriteVtuGeo(const std::string& datafile_name,
    const bool bool_binary, const bool bool_overlap,
    const DefAmrUint overlap_flag, const DefAmrIndexUint dims,
    const std::vector<DefAmrIndexUint>& vec_geo_in_one_vtu,
    const std::array<DefReal, 3>& grid_offset,
    OutputDataFormat& output_data_format,
    CriterionManager& criterion_manager) {
    FILE* fp = nullptr;
    std::string str_temp;
    errno_t err = fopen_s(&fp, (datafile_name + ".vtu").c_str(), "w");
    if (!fp) {
        int rank_id = 0;
#ifdef ENABLE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
        LogError("File on node " + std::to_string(rank_id)
            + " was not opened for writing vtu data in WriteVtuGeo.");
    } else {
        fprintf_s(fp, "<?xml version=\"1.0\"?>\n");
        str_temp.assign("<VTKFile type=\"UnstructuredGrid\""
            " version=\"0.1\" byte_order=\"" + str_vtk_byte_order_ + "\">\n");
        fprintf_s(fp, str_temp.c_str());
        fprintf_s(fp, " <UnstructuredGrid>\n");

        for (const auto& iter_geo : vec_geo_in_one_vtu) {
            WriteGeometryPieces(fp, bool_binary, dims, grid_offset,
                output_data_format,
                *criterion_manager.vec_ptr_geometries_.at(iter_geo));
        }

        fprintf_s(fp, " </UnstructuredGrid>\n");
        fprintf_s(fp, "</VTKFile>");
        fclose(fp);
    }
}
/**
* @brief   function to write xml formatted unstructed mesh (.vtu)
* @param[in]  program_name name of the program
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  bool_overlap write overlapping region or not
* @param[in]  overlap_flag flag indicates nodes in the overlapping region
* @param[in]  vec_level_in_one_vtu levels of grid write in this vtu file
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  grid_manager class to manage grid information
*/
void VtkWriterManager::WriteVtuGrid(const std::string& datafile_name,
    const bool bool_binary, const bool bool_overlap,
    const DefAmrUint overlap_flag,
    const std::vector<DefAmrIndexUint>& vec_level_in_one_vtu,
    OutputDataFormat& output_data_format,
    GridManagerInterface& grid_manager) {
    FILE* fp = nullptr;
    std::string str_temp;
    errno_t err = fopen_s(&fp, (datafile_name + ".vtu").c_str(), "w");
    if (!fp) {
        int rank_id = 0;
#ifdef ENABLE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
        LogError("File on node " + std::to_string(rank_id)
            + " was not opened for writing vtu data in WriteVtuGrid.");
    } else {
        fprintf_s(fp, "<?xml version=\"1.0\"?>\n");
        str_temp.assign("<VTKFile type=\"UnstructuredGrid\""
            " version=\"0.1\" byte_order=\"" + str_vtk_byte_order_ + "\">\n");
        fprintf_s(fp, str_temp.c_str());
        fprintf_s(fp, " <UnstructuredGrid>\n");

        for (const auto& iter_level : vec_level_in_one_vtu) {
            WriteGridPieces(fp, bool_binary, bool_overlap, overlap_flag,
                *grid_manager.vec_ptr_grid_info_.at(iter_level),
                output_data_format, grid_manager);
        }

        fprintf_s(fp, " </UnstructuredGrid>\n");
        fprintf_s(fp, "</VTKFile>");
        fclose(fp);
    }
}
/**
* @brief   function to write pieces of grids in vtu
* @param[in]  fp   pointer to output file
* @param[in]  vec_vtu_file_name name of vtu files
* @param[in]  output_data_format output data (real or integer) format
*/
void VtkWriterManager::WritePvtu(FILE* const fp,
    const std::vector<std::string> vec_vtu_file_name,
    OutputDataFormat& output_data_format) {
    std::string str_temp;
    fprintf_s(fp, "<?xml version=\"1.0\"?>\n");
    str_temp.assign("<VTKFile type=\"PUnstructuredGrid\""
        " version=\"0.1\" byte_order=\"" + str_vtk_byte_order_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, " <PUnstructuredGrid>\n");

    //// node data
    fprintf_s(fp, "    <PPointData>\n");
    // flag_satus
    str_temp.assign("      <PDataArray NumberOfComponents=\"1\" type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"flag_status\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "      </PDataArray>\n");
    // vtkGhostType
    str_temp.assign("      <PDataArray type=\"UInt8\" Name=\"vtkGhostType\""
        " format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "      </PDataArray>\n");
    fprintf_s(fp, "    </PPointData>\n");

    fprintf_s(fp, "   <PPoints>\n");
    // node coordinates
    str_temp.assign("    <PDataArray type=\""
        + output_data_format.output_real_.format_name_
        + "\" NumberOfComponents=\"3\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "    </PDataArray>\n");
    fprintf_s(fp, "   </PPoints>\n");

    fprintf_s(fp, "   <PCells>\n");
    // cell connectivity
    str_temp.assign("    <PDataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"connectivity\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "    </PDataArray>\n");
    // cell offset
    str_temp.assign("    <PDataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"offsets\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "    </PDataArray>\n");
    // cell types
    str_temp.assign("    <PDataArray type=\"UInt8\" Name=\"types\" "
        "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    fprintf_s(fp, "    </PDataArray>\n");

    fprintf_s(fp, "   </PCells>\n");

    for (const auto& iter : vec_vtu_file_name) {
        str_temp.assign("    <Piece Source=\""
            + iter + ".vtu" + "\"/>\n");
        fprintf_s(fp, str_temp.c_str());
    }

    fprintf_s(fp, " </PUnstructuredGrid>\n");
    fprintf_s(fp, " </VTKFile>\n");
}
/**
* @brief   function to write pieces of grids in vtu
* @param[in]  fp   pointer to output file
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  bool_overlap write overlapping region or not
* @param[in]  overlap_flag flag indicates nodes in the overlapping region
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  grid_manager class to manage grid information
*/
void VtkWriterManager::WriteGridPieces(FILE* const fp, const bool bool_binary,
    const bool bool_overlap, const DefAmrUint overlap_flag,
    const GridInfoInterface& grid_info,
    OutputDataFormat& output_data_format,
    GridManagerInterface& grid_manager) {

    std::string str_temp;
    DefAmrIndexUint dims = grid_manager.k0GridDims_;
    std::array<DefSizet, 2> node_and_cell_number;
    DefMap<DefSizet> map_node_index;
    if (dims == 2) {
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
        const GridManager2D& grid_manager2d =
            dynamic_cast<GridManager2D&> (grid_manager);
        node_and_cell_number = CalculateNumOfGridCells(
            bool_overlap, overlap_flag,
            grid_manager2d, grid_info.map_grid_node_, &map_node_index);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else {
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
        const GridManager3D& grid_manager3d =
            dynamic_cast<GridManager3D&> (grid_manager);
        node_and_cell_number = CalculateNumOfGridCells(
            bool_overlap, overlap_flag,
            grid_manager3d, grid_info.map_grid_node_, &map_node_index);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }
    str_temp.assign("  <Piece NumberOfPoints=\""
        + std::to_string(node_and_cell_number[0])
        + "\" NumberOfCells=\"" + std::to_string(node_and_cell_number[1])
        + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    fprintf_s(fp, "    <PointData>\n");
    // write node flag
    WriteGridNodeFlagStatus(fp, bool_binary, output_data_format,
        map_node_index, grid_info.map_grid_node_);
    WriteGridNodeVtkVisualization(fp, bool_binary, overlap_flag,
        output_data_format, map_node_index, grid_info.map_grid_node_);
    fprintf_s(fp, "    </PointData>\n");

    // write coordinates of vertices
    fprintf_s(fp, "   <Points>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_real_.format_name_
        + "\" NumberOfComponents=\"3\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    if (dims == 2) {
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
        const GridManager2D& grid_manager2d =
            dynamic_cast<GridManager2D&> (grid_manager);
        WriteGridCoordinates(fp, bool_binary,
            output_data_format, grid_manager2d, grid_info, &map_node_index);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else {
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
        const GridManager3D& grid_manager3d =
            dynamic_cast<GridManager3D&> (grid_manager);
        WriteGridCoordinates(fp, bool_binary,
            output_data_format, grid_manager3d, grid_info, &map_node_index);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }
    fprintf_s(fp, "    </DataArray>\n");
    fprintf_s(fp, "   </Points>\n");

    // write cell connectivity
    fprintf_s(fp, "   <Cells>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"connectivity\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    if (dims == 2) {
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
        const SFBitsetAux2D& bitset_aux2d =
            dynamic_cast<SFBitsetAux2D&> (grid_manager);
        WriteGridCellConnectivity(fp, bool_binary,
            output_data_format, bitset_aux2d, map_node_index);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else {
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
        const SFBitsetAux3D& bitset_aux3d =
            dynamic_cast<SFBitsetAux3D&> (grid_manager);
        WriteGridCellConnectivity(fp, bool_binary,
            output_data_format, bitset_aux3d, map_node_index);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }
    fprintf_s(fp, "    </DataArray>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"offsets\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    // write cell offsets
    WriteGridCellOffsets(fp, bool_binary, grid_manager.k0GridDims_,
        node_and_cell_number[1], output_data_format);
    fprintf_s(fp, "    </DataArray>\n");
    str_temp.assign("    <DataArray type=\"UInt8\" Name=\"types\" "
        "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    // write cell types
    WriteGridCellTypes(fp, bool_binary,
        grid_manager.k0GridDims_, node_and_cell_number[1]);
    fprintf_s(fp, "    </DataArray>\n");
    fprintf_s(fp, "   </Cells>\n");

    fprintf_s(fp, "  </Piece>\n");
}
/**
* @brief   function to write offset for each cell
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format.
* @param[in]  dims dimension of mesh.
* @param[in]  num_cell   number of cells.
* @param[in]  output_data_format output data (real or integer) format.
*/
void VtkWriterManager::WriteGridCellOffsets(FILE* const fp,
    const bool bool_binary, const DefAmrIndexUint dims,
    const DefSizet num_cell, OutputDataFormat& output_data_format) {

    DefAmrUint num_offsets = 4;
    if (dims == 2) {
        num_offsets = 4;  // VTK_PIXEL
    } else if (dims == 3) {
        num_offsets = 8;  // VTK_VOXEL
    }
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (DefSizet i = 1; i < num_cell + 1; ++i) {
            base64_instance_.AddToVecChar(
                output_data_format.output_uint_.CastType(
                    i * num_offsets), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "   "
            + output_data_format.output_uint_.printf_format_ + "\n";
        for (DefSizet i = 1; i < num_cell + 1; ++i) {
            fprintf_s(fp, str_format.c_str(), i * num_offsets);
        }
    }
}
/**
* @brief   function to write cell type.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  dims dimension
* @param[in]  num_cell   number of cells.
*/
void VtkWriterManager::WriteGridCellTypes(FILE* const fp,
    bool bool_binary, const DefAmrIndexUint dims, const DefSizet num_cell) {

    std::uint8_t cell_type = 8;
    if (dims == 2) {
        cell_type = 8;  // VTK_PIXEL
    } else if (dims == 3) {
        cell_type = 11;  // VTK_VOXEL
    }
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (DefSizet i = 0; i < num_cell; ++i) {
            base64_instance_.AddToVecChar(cell_type, &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (DefAmrIndexLUint i = 0; i < num_cell; ++i) {
            fprintf_s(fp, "   %u\n", cell_type);
        }
    }
}
/**
* @brief   function to write cell type.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  map_node_index    indices of nodes for unstructure grid.
* @param[in]  map_grid_node    grid nodes at the same refinement level.
*/
void VtkWriterManager::WriteGridNodeFlagStatus(FILE* const fp,
    const bool bool_binary, OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index,
    const DefMap<GridNode>& map_grid_node) {
    std::string str_temp;
    str_temp.assign("      <DataArray NumberOfComponents=\"1\" type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"flag_status\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            base64_instance_.AddToVecChar(
                output_data_format.output_uint_.CastType(
                    map_grid_node.at(iter->first).flag_status_), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "     "
            + output_data_format.output_uint_.printf_format_;
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            fprintf_s(fp, "  ");
            fprintf_s(fp, str_format.c_str(),
                map_grid_node.at(iter->first).flag_status_);
            fprintf_s(fp, "\n");
        }
    }
    fprintf_s(fp, "      </DataArray>\n");
}
/**
* @brief   function to write cell type.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format.
* @param[in]  flag_overlap_visual flag for vtk visualization.
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  map_node_index    indices of nodes for unstructure grid.
* @param[in]  map_grid_node    grid nodes at the same refinement level.
*/
void VtkWriterManager::WriteGridNodeVtkVisualization(
    FILE* const fp, const bool bool_binary,
    const DefAmrUint flag_overlap_visual,
    OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index,
    const DefMap<GridNode>& map_grid_node) {
    std::string str_temp;
    str_temp.assign("      <DataArray type=\"UInt8\" Name=\"vtkGhostType\" "
        "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
    uint8_t vtk_normal = kVtkNormalPoint_, vtk_overlap = kVtkDuplicatedPoint_;
    fprintf_s(fp, str_temp.c_str());
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            if ((map_grid_node.at(iter->first).flag_status_
                & flag_overlap_visual) == 0) {
                base64_instance_.AddToVecChar(vtk_normal, &vec_uint8);
            } else {
                base64_instance_.AddToVecChar(vtk_overlap, &vec_uint8);
            }
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            if ((map_grid_node.at(iter->first).flag_status_
                & flag_overlap_visual) == 0) {
                fprintf_s(fp, "   %u\n", vtk_normal);
            } else {
                fprintf_s(fp, "   %u\n", vtk_overlap);
            }
        }
    }
    fprintf_s(fp, "      </DataArray>\n");
}
/**
* @brief   function to write pieces of geometries in vtu.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  dims dimension of the mesh
* @param[in]  grid_offset offset coordinates of the mesh.
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  criterion_manager class manager criterion
*/
void VtkWriterManager::WriteGeometryPieces(
    FILE* const fp, const bool bool_binary,
    const DefAmrIndexUint dims, const std::array<DefReal, 3>& grid_offset,
    OutputDataFormat& output_data_format,
    GeometryInfoInterface& geo_info) {
    std::string str_temp;
    DefSizet cell_number = CalculateNumOfGeometryCells(geo_info);
    str_temp.assign("  <Piece NumberOfPoints=\""
        + std::to_string(geo_info.GetNumOfGeometryPoints())
        + "\" NumberOfCells=\"" + std::to_string(cell_number)
        + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    // write vertex data
    fprintf_s(fp, "   <Points>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_real_.format_name_
        + "\" NumberOfComponents=\"3\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    DefSizet num_points = 0;
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
    if (dims == 2) {
        const GeometryInfo2DInterface& geo_info_2d =
            dynamic_cast<GeometryInfo2DInterface&> (geo_info);
        num_points = WriteGeometryCoordinates(fp, bool_binary, dims,
            grid_offset, output_data_format, geo_info_2d);
    }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
    if (dims == 3) {
        const GeometryInfo3DInterface& geo_info_3d =
            dynamic_cast<GeometryInfo3DInterface&> (geo_info);
        num_points = WriteGeometryCoordinates(fp, bool_binary, dims,
            grid_offset, output_data_format, geo_info_3d);
    }
#endif  // DEBUG_DISABLE_3D_FUNCTIONS

    fprintf_s(fp, "    </DataArray>\n");
    fprintf_s(fp, "   </Points>\n");
    fprintf_s(fp, "   <Cells>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"connectivity\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());
    // write cell connectiveity
    DefSizet num_cell = 0;
    switch (geo_info.geometry_cell_type_) {
    case   EGeometryCellType::kPolyLine: {
        WriteGeometryCellConnectivitykPolyLine(
            fp, bool_binary, num_points, output_data_format);
        num_cell = num_points - 1;
        break;
    }
    default:
        break;
    }

    fprintf_s(fp, "    </DataArray>\n");
    str_temp.assign("    <DataArray type=\""
        + output_data_format.output_uint_.format_name_
        + "\" Name=\"offsets\" format=\""
        + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    // write cell offsets
    WriteGeometryCellOffset(fp, bool_binary, num_cell,
        geo_info.geometry_cell_type_, output_data_format);

    fprintf_s(fp, "    </DataArray>\n");
    str_temp.assign("    <DataArray type=\"UInt8\" Name=\"types\" "
        "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    // write cell types
    WriteGeometryCellType(
        fp, bool_binary, num_points, geo_info.geometry_cell_type_);

    fprintf_s(fp, "    </DataArray>\n");
    fprintf_s(fp, "   </Cells>\n");
    fprintf_s(fp, "  </Piece>\n");
}
/**
* @brief   function to calculate number of cells.
*          according to geometry_cell_type_.
* @param[in]  geo_info   instance of geometry information.
* @return  number of cells.
*/
DefSizet VtkWriterManager::CalculateNumOfGeometryCells(
    const GeometryInfoInterface& geo_info) const {
    DefSizet num_cell = 0;
    switch (geo_info.geometry_cell_type_) {
    case   EGeometryCellType::kPolyLine:
        num_cell = geo_info.GetNumOfGeometryPoints() - 1;
        break;
    default:
        LogError("GeometryOutputType is undefined for calculing cell"
            " numbers in CalculateGemetryCells.");
        break;
    }
    return num_cell;
}
/**
* @brief   function to write connection informtion for polyline
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  num_points number of geometry points
* @param[in]  output_data_format output data (real or integer) format
*/
void VtkWriterManager::WriteGeometryCellConnectivitykPolyLine(FILE* const fp,
    const bool bool_binary, const DefSizet num_points,
    OutputDataFormat& output_data_format) {
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (DefSizet i_vertex = 0; i_vertex < num_points - 1; ++i_vertex) {
            base64_instance_.AddToVecChar(
                output_data_format.output_uint_.CastType(
                    i_vertex), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_uint_.CastType(
                    i_vertex + 1), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "   "
            + output_data_format.output_uint_.printf_format_ + " "
            + output_data_format.output_uint_.printf_format_ + "\n";
        for (DefSizet i_vertex = 0; i_vertex < num_points - 1; ++i_vertex) {
            fprintf_s(fp, str_format.c_str(), i_vertex,
                i_vertex + 1);
        }
    }
}
/**
* @brief   function to write offset of each cell for polyline.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format.
* @param[in]  geometry_cell_type cell type of geometry.
* @param[in]  num_cell number of geometry cells.
* @param[in]  output_data_format output data (real or integer) format.
*/
void VtkWriterManager::WriteGeometryCellOffset(FILE* const fp,
    const bool bool_binary, const DefSizet num_cell,
    const EGeometryCellType geometry_cell_type,
    OutputDataFormat& output_data_format) {
    DefSizet num_offset = 0;
    switch (geometry_cell_type) {
    case EGeometryCellType::kPolyLine:
        num_offset = 2;
        break;
    case EGeometryCellType::kTriangle:
        num_offset = 3;
        break;
    default:
        break;
    }

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (DefSizet i_vertex = 1; i_vertex < num_cell + 1; ++i_vertex) {
            base64_instance_.AddToVecChar(
                output_data_format.output_uint_.CastType(
                    i_vertex * num_offset), &vec_uint8);;
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "   "
            + output_data_format.output_uint_.printf_format_ + "\n";
        for (DefSizet i_vertex = 1; i_vertex < num_cell + 1; ++i_vertex) {
            fprintf_s(fp, str_format.c_str(), i_vertex * num_offset);
        }
    }
}
/**
* @brief   function to write cell type.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  num_points number of geometry points
* @param[in]  geometry_cell_type cell type of geometry
*/
void VtkWriterManager::WriteGeometryCellType(
    FILE* const fp, const bool bool_binary, const DefSizet num_points,
    const EGeometryCellType geometry_cell_type) {

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (DefSizet i_vertex = 0; i_vertex < num_points - 1; ++i_vertex) {
            base64_instance_.AddToVecChar(static_cast<std::uint8_t>(
                geometry_cell_type), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (DefSizet i_vertex = 0; i_vertex < num_points - 1; ++i_vertex) {
            fprintf_s(fp, "   %u\n",
                static_cast<std::uint8_t>(geometry_cell_type));
        }
    }
}
#ifndef DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief   function to calculate number of cells
*          according to geometry_cell_type_.
* @param[in]  sfbitset_aux2d   class to manage functions of spacing
*                          filling code related manipulations.
* @param[in]  map_grid_node   grid nodes at the same refinement level.
* @param[out]  ptr_map_node_index   indices of grid nodes
*                           in writing unstructure grid.
* @return  number of nodes and number of cells.
*/
std::array<DefSizet, 2> VtkWriterManager::CalculateNumOfGridCells(
    const bool bool_overlap, const DefAmrUint overlap_flag,
    const SFBitsetAux2D& sfbitset_aux2d,
    const DefMap<GridNode>& map_grid_node,
    DefMap<DefSizet>* ptr_map_node_index) const {

    DefSizet num_cell = 0;
    std::array<DefSFBitset, 4> bitset_cell;
    bool no_overlap;
    for (auto iter = map_grid_node.begin();
        iter != map_grid_node.end(); ++iter) {
        if (bool_overlap ||  // points not in overlapping regions
            (iter->second.flag_status_ & overlap_flag) == 0) {
            // record the node indices in map_grid_node_
            ptr_map_node_index->insert({ iter->first, 0 });
            // check if the node is at the lower corner of a cell
            if (sfbitset_aux2d.SFBitsetBelongToOneCell(
                iter->first, map_grid_node, &bitset_cell)) {
                no_overlap = true;
                for (const auto& iter_cell : bitset_cell) {
                    if (!bool_overlap
                        && (map_grid_node.at(iter_cell).flag_status_
                        & overlap_flag) != 0) {
                        no_overlap = false;
                        break;
                    }
                }
                if (no_overlap) {
                    ++num_cell;
                }
            }
        }
    }
    DefSizet num_nodes = ptr_map_node_index->size();
    return { num_nodes, num_cell };
}
/**
* @brief   function to write grid coordinates.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  grid_manager2d class to manage 3d grid information
* @param[in]  grid_info   instance of grid information.
* @param[in]  map_node_index    indices of nodes for unstructure grid.
*/
void VtkWriterManager::WriteGridCoordinates(FILE* const fp,
    const bool bool_binary,
    OutputDataFormat& output_data_format,
    const GridManager2D& grid_manager2d,
    const GridInfoInterface& grid_info,
    DefMap<DefSizet>* const ptr_map_node_index) {
    std::string str_format = "   "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + "\n";
    DefReal x_offset = 0., y_offset = 0., z_offset = 0.;
    std::array<DefReal, 2> grid_space_2d, coordi_2d;

    DefSizet i_nodes = 0;
    x_offset = grid_manager2d.k0RealOffset_.at(kXIndex);
    y_offset = grid_manager2d.k0RealOffset_.at(kYIndex);
    grid_space_2d = { grid_info.grid_space_[kXIndex],
        grid_info.grid_space_[kYIndex] };

    std::vector<DefReal> coordi(3, 0.);
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = ptr_map_node_index->begin();
            iter != ptr_map_node_index->end(); ++iter) {
            iter->second = i_nodes;
            // compute coordinates
            grid_manager2d.SFBitsetComputeCoordinate(iter->first,
                grid_space_2d, &coordi_2d);
            coordi = { coordi_2d[kXIndex], coordi_2d[kYIndex], 0 };
            // convert to uint8
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kXIndex] - x_offset), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kYIndex] - y_offset), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kZIndex] - z_offset), &vec_uint8);

            ++i_nodes;
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (auto iter = ptr_map_node_index->begin();
            iter != ptr_map_node_index->end(); ++iter) {
            iter->second = i_nodes;
            grid_manager2d.SFBitsetComputeCoordinate(iter->first,
                grid_space_2d, &coordi_2d);
            coordi = { coordi_2d[kXIndex], coordi_2d[kYIndex], 0 };
            fprintf_s(fp, str_format.c_str(),
                coordi[kXIndex] - x_offset,
                coordi[kYIndex] - y_offset,
                coordi[kZIndex] - z_offset);

            ++i_nodes;
        }
    }
}
/**
* @brief   function to write connectivity relation of cells
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  bitset_aux2d   class to manage functions of spacing
*                          filling code related manipulations.
* @param[in]  map_grid_node   grid nodes
* @param[in]  map_node_index   indices of nodes for unstructure grid.
*/
void VtkWriterManager::WriteGridCellConnectivity(FILE* const fp,
    const bool bool_binary,
    OutputDataFormat& output_data_format,
    const SFBitsetAux2D& bitset_aux2d,
    const DefMap<DefSizet>& map_node_index) {
    DefSizet size_vec = TwoPowerN(static_cast<DefSizet>(2));
    std::array<DefSFBitset, 4> bitset_cell_2d;
    // write connection
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            // check if the node is at the lower corner of a cell
            if (bitset_aux2d.SFBitsetBelongToOneCell(
                iter->first, map_node_index, &bitset_cell_2d)) {
                for (const auto& iter_sfbitset : bitset_cell_2d) {
                    base64_instance_.AddToVecChar(
                        output_data_format.output_uint_.CastType(
                            map_node_index.at(iter_sfbitset)), &vec_uint8);
                }
            }
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = " "
            + output_data_format.output_uint_.printf_format_;
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            if (bitset_aux2d.SFBitsetBelongToOneCell(
                iter->first, map_node_index, &bitset_cell_2d)) {
                fprintf_s(fp, "  ");
                for (const auto& iter_sfbitset : bitset_cell_2d) {
                    fprintf_s(fp, str_format.c_str(),
                        map_node_index.at(iter_sfbitset));
                }
                fprintf_s(fp, "\n");
            }
        }
    }
}
/**
* @brief   function to write geometry coordinates.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  dims dimension
* @param[in]  grid_offset offset coordinates of the mesh
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  geo_info class to manage geometry information
*/
DefSizet VtkWriterManager::WriteGeometryCoordinates(FILE* const fp,
    const bool bool_binary, const DefAmrIndexUint dims,
    const std::array<DefReal, 3>& grid_offset,
    OutputDataFormat& output_data_format,
    const GeometryInfo2DInterface& geo_info) {

    std::string str_format = "   "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + "\n";
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (const auto& iter : geo_info.coordinate_origin_) {
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    iter.coordinate.at(kXIndex)
                    - grid_offset.at(kXIndex)), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    iter.coordinate.at(kYIndex)
                    - grid_offset.at(kYIndex)), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(0.), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (const auto& iter : geo_info.coordinate_origin_) {
            fprintf_s(fp, str_format.c_str(),
                iter.coordinate.at(kXIndex) - grid_offset.at(kXIndex),
                iter.coordinate.at(kYIndex) - grid_offset.at(kYIndex),
                0.);
        }
    }
    return geo_info.coordinate_origin_.size();
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to calculate number of cells
*          according to geometry_cell_type_.
* @param[in]  map_grid_node   grid nodes at the same refinement level.
* @param[in]  bool_overlap write overlapping region or not
* @param[in]  overlap_flag flag indicates nodes in the overlapping region
* @param[out]  ptr_map_node_index   indices of grid nodes
*                           for unstructure grid.
* @return  number of nodes and number of cells.
*/
std::array<DefSizet, 2> VtkWriterManager::CalculateNumOfGridCells(
    const bool bool_overlap, const DefAmrUint overlap_flag,
    const SFBitsetAux3D& sfbitset_aux3d,
    const DefMap<GridNode>& map_grid_node,
    DefMap<DefSizet>* ptr_map_node_index) const {

    DefAmrIndexLUint num_cell = 0;
    DefAmrIndexLUint node_index = 0;
    std::array<DefSFBitset, 8> bitset_cell;

    bool no_overlap;
    for (auto iter = map_grid_node.begin();
        iter != map_grid_node.end(); ++iter) {
        if (bool_overlap ||  // points not in overlapping regions
            (iter->second.flag_status_ & overlap_flag) == 0) {
            // record the node indices in map_grid_node_
            ptr_map_node_index->insert({ iter->first, 0 });
            // check if the node is at the lower corner of a cell
            if (sfbitset_aux3d.SFBitsetBelongToOneCell(
                iter->first, map_grid_node, &bitset_cell)) {
                no_overlap = true;
                for (const auto& iter_cell : bitset_cell) {
                    if (!bool_overlap
                        && (map_grid_node.at(iter_cell).flag_status_
                            & overlap_flag) != 0) {
                        no_overlap = false;
                        break;
                    }
                }
                if (no_overlap) {
                    ++num_cell;
                }
            }
        }
    }
    DefSizet num_nodes = ptr_map_node_index->size();
    return { num_nodes, num_cell };
}
/**
* @brief   function to write grid coordinates.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  grid_manager2d class to manage 3d grid information
* @param[in]  grid_info   instance of grid information.
* @param[in]  map_node_index    indices of nodes for unstructure grid.
*/
void VtkWriterManager::WriteGridCoordinates(FILE* const fp,
    const bool bool_binary,
    OutputDataFormat& output_data_format,
    const GridManager3D& grid_manager3d,
    const GridInfoInterface& grid_info,
    DefMap<DefSizet>* const ptr_map_node_index) {
    std::string str_format = "   "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + "\n";
    DefReal x_offset = 0., y_offset = 0., z_offset = 0.;
    std::array<DefReal, 3> grid_space_3d, coordi_3d;
    x_offset = grid_manager3d.k0RealOffset_.at(kXIndex);
    y_offset = grid_manager3d.k0RealOffset_.at(kYIndex);
    z_offset = grid_manager3d.k0RealOffset_.at(kZIndex);
    grid_space_3d = { grid_info.grid_space_[kXIndex],
        grid_info.grid_space_[kYIndex],  grid_info.grid_space_[kZIndex] };
    std::vector<DefReal> coordi(3, 0.);
    DefSizet i_node = 0;
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = ptr_map_node_index->begin();
            iter != ptr_map_node_index->end(); ++iter) {
            iter->second = i_node;
            // compute coordinates
                grid_manager3d.SFBitsetComputeCoordinate(iter->first,
                    grid_space_3d, &coordi_3d);
                coordi = { coordi_3d[kXIndex], coordi_3d[kYIndex],
                    coordi_3d[kZIndex] };
            // convert to uint8
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kXIndex] - x_offset), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kYIndex] - y_offset), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    coordi[kZIndex] - z_offset), &vec_uint8);

            ++i_node;
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (auto iter = ptr_map_node_index->begin();
            iter != ptr_map_node_index->end(); ++iter) {
            iter->second = i_node;
                grid_manager3d.SFBitsetComputeCoordinate(iter->first,
                    grid_space_3d, &coordi_3d);
                coordi = { coordi_3d[kXIndex], coordi_3d[kYIndex],
                    coordi_3d[kZIndex] };
            fprintf_s(fp, str_format.c_str(),
                coordi[kXIndex] - x_offset,
                coordi[kYIndex] - y_offset,
                coordi[kZIndex] - z_offset);

            ++i_node;
        }
    }
}
/**
* @brief   function to write connectivity relation of cells
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  bitset_aux3d   class to manage functions of spacing
*                          filling code related manipulations.
* @param[in]  map_node_index   indices of nodes for unstructure grid.
*/
void VtkWriterManager::WriteGridCellConnectivity(FILE* const fp,
    const bool bool_binary,
    OutputDataFormat& output_data_format,
    const SFBitsetAux3D& bitset_aux3d,
    const DefMap<DefSizet>& map_node_index) {
    DefSizet size_vec = TwoPowerN(static_cast<DefSizet>(3));
    std::array<DefSFBitset, 8> bitset_cell_3d;

    // write connection
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            // check if the node is at the lower corner of a cell
            if (bitset_aux3d.SFBitsetBelongToOneCell(
                iter->first, map_node_index, &bitset_cell_3d)) {
                for (const auto& iter_sfbitset : bitset_cell_3d) {
                    base64_instance_.AddToVecChar(
                        output_data_format.output_uint_.CastType(
                            map_node_index.at(iter_sfbitset)), &vec_uint8);
                }
            }
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = " "
            + output_data_format.output_uint_.printf_format_;
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            // check if the node is at the lower corner of a cell
            if (bitset_aux3d.SFBitsetBelongToOneCell(
                iter->first, map_node_index, &bitset_cell_3d)) {
                fprintf_s(fp, "  ");
                for (const auto& iter_sfbitset : bitset_cell_3d) {
                    fprintf_s(fp, str_format.c_str(),
                        map_node_index.at(iter_sfbitset));
                }
                fprintf_s(fp, "\n");
            }
        }
    }
}
/**
* @brief   function to write geometry coordinates.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in bindary or ascii format
* @param[in]  dims dimension
* @param[in]  grid_offset offset coordinates of the mesh
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  geo_info   instance of geometry information.
*/
DefSizet VtkWriterManager::WriteGeometryCoordinates(FILE* const fp,
    const bool bool_binary, const DefAmrIndexUint dims,
    const std::array<DefReal, 3>& grid_offset,
    OutputDataFormat& output_data_format,
    const GeometryInfo3DInterface& geo_info) {

    std::string str_format = "   "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + " "
        + output_data_format.output_real_.printf_format_ + "\n";
    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (const auto& iter : geo_info.coordinate_origin_) {
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    iter.coordinate.at(kXIndex)
                    - grid_offset.at(kXIndex)), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    iter.coordinate.at(kYIndex)
                    - grid_offset.at(kYIndex)), &vec_uint8);
            base64_instance_.AddToVecChar(
                output_data_format.output_real_.CastType(
                    iter.coordinate.at(kZIndex)
                    - grid_offset.at(kZIndex)), &vec_uint8);
        }
        base64_instance_.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        for (const auto& iter : geo_info.coordinate_origin_) {
            fprintf_s(fp, str_format.c_str(),
                iter.coordinate.at(kXIndex) - grid_offset.at(kXIndex),
                iter.coordinate.at(kYIndex) - grid_offset.at(kYIndex),
                iter.coordinate.at(kZIndex) - grid_offset.at(kZIndex));
        }
    }
    return geo_info.coordinate_origin_.size();
}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
