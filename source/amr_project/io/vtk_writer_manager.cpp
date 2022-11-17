////  Copyright (c) 2022, Zhengliang Liu
////  All rights reserved
//
///**
//* @file vtk_writer_manager.cpp
//* @author Zhengliang Liu
//* @brief functions used to write vtk data.
//* @date  2022-8-04
//*/
//#include <typeinfo>
//#include "auxiliary_inline_func.h"
//#include "io/log_write.h"
//#include "io/io_manager.h"
//#include "criterion/criterion_manager.h"
//#include "grid/grid_manager.h"
//#include "mpi/mpi_manager.h"
//namespace rootproject {
//namespace amrproject {
//namespace io {
///**
//* @brief  funciton to convert data to base64
//* @param[in]  char in   vertexer vertex to output file
//*/
//void Base64Utility::Encode(std::vector<uint8_t>* ptr_vec_uint8,
//    std::vector<uint8_t>* ptr_vec_base64) {
//    // number of bytes of the output data
//    const uint32_t size_vec = static_cast<uint32_t>(ptr_vec_uint8->size());
//    ptr_vec_uint8->insert(ptr_vec_uint8->begin(),
//        reinterpret_cast<const uint8_t*> (&size_vec),
//        reinterpret_cast<const uint8_t*> (&size_vec) + sizeof(uint32_t));
//    DefSizet uint8_length = ptr_vec_uint8->size();
//    DefSizet base64_length = uint8_length / 3 * 4;
//    if (base64_length % 3 > 0) {
//        base64_length += 4;
//    }
//    ptr_vec_base64->resize(base64_length);
//
//    DefSizet i, index = 0;
//    for (i = 0; i + 3 <= uint8_length; i += 3) {
//        ptr_vec_base64->at(index++) =
//            kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
//        ptr_vec_base64->at(index++) =
//            kEncodeTable_[((ptr_vec_uint8->at(i) << 4) & 0x30)
//            | ((ptr_vec_uint8->at(i + 1) >> 4) & 0x0F)];
//        ptr_vec_base64->at(index++) =
//            kEncodeTable_[((ptr_vec_uint8->at(i + 1) << 2) & 0x3C)
//            | (ptr_vec_uint8->at(i + 2) >> 6) & 0x03];
//        ptr_vec_base64->at(index++) =
//            kEncodeTable_[ptr_vec_uint8->at(i + 2) & 0x3F];
//    }
//
//    if (i < uint8_length) {
//        DefSizet bytes_remain = uint8_length - i;
//        if (bytes_remain == 1) {  // single byte
//            ptr_vec_base64->at(index++) =
//                kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
//            ptr_vec_base64->at(index++) =
//                kEncodeTable_[(ptr_vec_uint8->at(i) << 4) & 0x30];
//            ptr_vec_base64->at(index++) = '=';
//            ptr_vec_base64->at(index++) = '=';
//        } else  {  // pair bytes
//            ptr_vec_base64->at(index++) =
//                kEncodeTable_[(ptr_vec_uint8->at(i) >> 2) & 0x3F];
//            ptr_vec_base64->at(index++) =
//                kEncodeTable_[((ptr_vec_uint8->at(i) << 4) & 0x30)
//                | ((ptr_vec_uint8->at(i + 1) >> 4)) & 0x0F];
//            ptr_vec_base64->at(index++) =
//                kEncodeTable_[(ptr_vec_uint8->at(i + 1) << 2) & 0x3C];
//            ptr_vec_base64->at(index++) = '=';
//        }
//    }
//}
///**
//* @brief  funciton to initialize vtk related options
//*/
//void VtkWriterManager::OptionInitial(const bool bool_binary) {
//    if (CheckIfLittleEndian()) {
//        str_vtk_byte_order_ = "LittleEndian";
//    } else {
//        str_vtk_byte_order_ = "BigEndian";
//    }
//    if (bool_binary) {
//        k0StrVtkAsciiOrBinary_ = "binary";
//    } else {
//        k0StrVtkAsciiOrBinary_ = "ascii";
//    }
//}
///**
//* @brief   function to write xml formatted unstructed mesh (.vtu)
//*/
//void VtkWriterManager::WriteVtu(const bool bool_binary,
//    OutputDataFormat& output_data_format,
//    std::shared_ptr<grid::GridManager> ptr_grid_manager,
//    std::shared_ptr<criterion::CriterionManager> ptr_criterion_manager) {
//    std::string datafile_name = "test";
//    FILE* fp = nullptr;
//    std::string str_temp;
//    errno_t err = fopen_s(&fp,  (datafile_name + ".vtu").c_str(), "w");
//    if (!fp) {
//        int rank_id = 0;
//#ifdef ENABLE_MPI
//        rank_id = mpi::MpiManager::GetInstance()->rank_id_;
//#endif  // ENABLE_MPI
//        io::LogError("File was not opened for writing vtu data in WriteVtu.");
//    } else {
//        fprintf_s(fp, "<?xml version=\"1.0\"?>\n");
//        str_temp.assign("<VTKFile type=\"UnstructuredGrid\""
//            " version=\"0.1\" byte_order=\"" + str_vtk_byte_order_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//        fprintf_s(fp, " <UnstructuredGrid>\n");
//
//        WriteGrid(fp, bool_binary, output_data_format, ptr_grid_manager);
//
//        WriteGeometry(fp, bool_binary, output_data_format, ptr_criterion_manager);
//
//        fprintf_s(fp, " </UnstructuredGrid>\n");
//        fprintf_s(fp, "</VTKFile>");
//        fclose(fp);
//    }
//}
///**
//* @brief   function to write pieces of grids in vtu
//* @param[in]  fp   vertexer vertex to output file
//*/
//void VtkWriterManager::WriteGrid(FILE* fp, const bool bool_binary,
//    OutputDataFormat& output_data_format,
//    std::shared_ptr<grid::GridManager> ptr_grid_manager) {
//    std::string str_temp;
//    std::array<DefLUint, 2> node_and_cell_number;
//    for (DefSizet i_level = 0;
//        i_level <= ptr_grid_manager->k0MaxLevel_; ++i_level) {
//        DefMap<DefLUint> map_node_index;
//        node_and_cell_number = CalculateNumOfGridCells(
//            ptr_grid_manager->k0GridDims_,
//            ptr_grid_manager->vec_ptr_grid_info_.at(i_level)->map_grid_node_,
//            &map_node_index);
//        str_temp.assign("  <Piece NumberOfvertexs=\""
//            + std::to_string(node_and_cell_number[0])
//            + "\" NumberOfCells=\"" + std::to_string(node_and_cell_number[1])
//            + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//        fprintf_s(fp, "   <vertexs>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_real_.format_name_
//            + "\" NumberOfComponents=\"3\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write vertex data
//        WirteGridCoordinates(fp, bool_binary,
//            output_data_format, ptr_grid_manager,
//            *(ptr_grid_manager->vec_ptr_grid_info_.at(i_level)));
//
//        fprintf_s(fp, "    </DataArray>\n");
//        fprintf_s(fp, "   </vertexs>\n");
//        fprintf_s(fp, "   <Cells>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_uint_.format_name_
//            + "\" Name=\"connectivity\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell connectiveity
//        WirteGridCellConnectivity(fp, bool_binary,
//            ptr_grid_manager->k0GridDims_, output_data_format,
//            ptr_grid_manager->vec_ptr_grid_info_.at(i_level)->map_grid_node_,
//            map_node_index);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_uint_.format_name_
//            + "\" Name=\"offsets\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell offsets
//        WirteGridCellOffsets(fp, bool_binary, ptr_grid_manager->k0GridDims_,
//            node_and_cell_number[1], output_data_format);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        str_temp.assign("    <DataArray type=\"UInt8\" Name=\"types\" "
//            "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell types
//        WirteGridCellTypes(fp, bool_binary,
//            ptr_grid_manager->k0GridDims_, node_and_cell_number[1]);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        fprintf_s(fp, "   </Cells>\n");
//        fprintf_s(fp, "  </Piece>\n");
//    }
//}
///**
//* @brief   function to calculate number of cells
//*          according to geometry_cell_type_.
//* @param[in]  map_grid_node   grid nodes at the same refinement level.
//* @param[out]  ptr_map_node_index   indices of grid nodes
//*                           in writing unstructure grid.
//* @return  number of nodes and number of cells.
//*/
//std::array<DefLUint, 2> VtkWriterManager::CalculateNumOfGridCells(
//    const DefUint dims,
//    const DefMap<grid::GridNode>& map_grid_node,
//    DefMap<DefLUint>* ptr_map_node_index) {
//
//    DefLUint num_cell = 0;
//    DefLUint node_index = 0;
//    std::array<DefSFBitset, 4> sfbitset2D(4, 0);
//
//    // vertexer to dimesional related function
//    bool(*ptrfunc_bitset_belong_to_one_cell)(
//        const DefSFBitset&, const DefMap<grid::GridNode>&,
//        std::vector<DefSFBitset>*) = nullptr;
//    if (dims == 2) {
//#ifndef DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            &grid::SFBitsetBelongToOneCell2D<grid::GridNode>;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    } else {
//#ifndef DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            &grid::SFBitsetBelongToOneCell3D<grid::GridNode>;
//    }
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//
//    for (auto iter = map_grid_node.begin();
//        iter != map_grid_node.end(); ++iter) {
//        // record the node indices in map_grid_node_
//        ptr_map_node_index->insert({ iter->first, node_index });
//        ++node_index;
//
//        // check if the node is at the lower corner of a cell
//        if (dims == 2) {
//            if (grid::SFBitsetBelongToOneCell2D<grid::GridNode>(
//                iter->first, map_grid_node, sfbitset2D)) {
//
//            }
//        }
//        if ((*ptrfunc_bitset_belong_to_one_cell)(
//            iter->first, map_grid_node, &vec_sfbitset)) {
//            ++num_cell;
//        }
//    }
//
//    return { node_index, num_cell };
//}
///**
//* @brief   function to write geometry coordinates.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information.
//*/
//void VtkWriterManager::WirteGridCoordinates(FILE* fp,
//    const bool bool_binary, OutputDataFormat& output_data_format,
//    std::shared_ptr<grid::GridManager> ptr_grid_manager,
//    const grid::GridInfoInterface& grid_info) {
//    std::string str_format = "   "
//        + output_data_format.output_real_.printf_format_ + " "
//        + output_data_format.output_real_.printf_format_ + " "
//        + output_data_format.output_real_.printf_format_ + "\n";
//    DefReal x_offset = ptr_grid_manager->k0RealOffest_.at(kXIndex),
//        y_offset = ptr_grid_manager->k0RealOffest_.at(kYIndex);
//    DefReal z_offset = 0.;
//    std::vector<DefReal> coordi(3, 0.);
//    if (ptr_grid_manager->k0GridDims_ == 3) {
//        z_offset = ptr_grid_manager->k0RealOffest_.at(kZIndex);
//    }
//
//    // vertexer to dimension specified function
//    void (*ptrfunc_sfbitset_compute_coordinate_)(
//        const DefSFBitset & bitset_in,
//        const std::vector<DefReal>&vec_grid_space,
//        std::vector<DefReal>*ptr_vec_coordi) = nullptr;
//    if (ptr_grid_manager->k0GridDims_ == 2) {
//#ifndef DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_sfbitset_compute_coordinate_ =
//            &grid::SFBitsetComputeCoordinate2D;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    } else {
//#ifndef DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_sfbitset_compute_coordinate_ =
//            &grid::SFBitsetComputeCoordinate3D;
//    }
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//
//    if (bool_binary) {
//        std::vector<uint8_t> vec_uint8{}, vec_base64{};
//        for (auto iter = grid_info.map_grid_node_.begin();
//            iter != grid_info.map_grid_node_.end(); ++iter) {
//            // compute coordinates
//            ptrfunc_sfbitset_compute_coordinate_(iter->first,
//                grid_info.grid_space_, &coordi);
//            // convert to uint8
//            base64_instance_.add_to_vec_char(
//                output_data_format.output_real_.cast_type(
//                    coordi[kXIndex] - x_offset), &vec_uint8);
//            base64_instance_.add_to_vec_char(
//                output_data_format.output_real_.cast_type(
//                    coordi[kYIndex] - y_offset), &vec_uint8);
//            base64_instance_.add_to_vec_char(
//                output_data_format.output_real_.cast_type(
//                    coordi[kZIndex] - z_offset), &vec_uint8);
//        }
//        base64_instance_.Encode(&vec_uint8, &vec_base64);
//        for (const auto& iter : vec_base64) {
//            fprintf_s(fp, "%c", iter);
//        }
//        fprintf_s(fp, "\n");
//    } else {
//        for (auto iter = grid_info.map_grid_node_.begin();
//            iter != grid_info.map_grid_node_.end(); ++iter) {
//            ptrfunc_sfbitset_compute_coordinate_(iter->first,
//                grid_info.grid_space_, &coordi);
//            fprintf_s(fp, str_format.c_str(),
//                coordi[kXIndex] - x_offset,
//                coordi[kYIndex] - y_offset,
//                coordi[kZIndex] - z_offset);
//        }
//    }
//}
///**
//* @brief   function to write connectivity relation of cells
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information
//*/
//void VtkWriterManager::WirteGridCellConnectivity(FILE* fp,
//    const bool bool_binary, const DefUint dims,
//    OutputDataFormat& output_data_format,
//    const DefMap<grid::GridNode>& map_grid_node,
//    const DefMap<DefLUint>& map_node_index) {
//
//    DefSizet size_vec = TwoPowerN(static_cast<DefSizet>(dims));
//    std::vector<DefSFBitset> vec_sfbitset(size_vec, 0);
//
//    // vertexer to dimesional related function
//    bool(*ptrfunc_bitset_belong_to_one_cell)(
//        const DefSFBitset&, const DefMap<grid::GridNode>&,
//        std::vector<DefSFBitset>*) = nullptr;
//    if (dims == 2) {
//#ifndef DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            &grid::SFBitsetBelongToOneCell2D<grid::GridNode>;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    } else {
//#ifndef DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            &grid::SFBitsetBelongToOneCell3D<grid::GridNode>;
//    }
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//
//    // write connection
//    if (bool_binary) {
//        std::vector<uint8_t> vec_uint8{}, vec_base64{};
//        for (auto iter = map_grid_node.begin();
//            iter != map_grid_node.end(); ++iter) {
//            // check if the node is at the lower corner of a cell
//            if ((ptrfunc_bitset_belong_to_one_cell)(
//                iter->first, map_grid_node, &vec_sfbitset)) {
//                for (const auto& iter_sfbitset : vec_sfbitset) {
//                    base64_instance_.add_to_vec_char(
//                        output_data_format.output_uint_.cast_type(
//                            map_node_index.at(iter_sfbitset)), &vec_uint8);
//                }
//            }
//        }
//        base64_instance_.Encode(&vec_uint8, &vec_base64);
//        for (const auto& iter : vec_base64) {
//            fprintf_s(fp, "%c", iter);
//        }
//        fprintf_s(fp, "\n");
//    } else {
//        std::string str_format = " "
//            + output_data_format.output_uint_.printf_format_;
//        for (auto iter = map_grid_node.begin();
//            iter != map_grid_node.end(); ++iter) {
//            // check if the node is at the lower corner of a cell
//            if ((ptrfunc_bitset_belong_to_one_cell)(
//                iter->first, map_grid_node, &vec_sfbitset)) {
//                fprintf_s(fp, "  ");
//                for (const auto& iter_sfbitset : vec_sfbitset) {
//                    fprintf_s(fp, str_format.c_str());
//                }
//                fprintf_s(fp, "\n");
//            }
//        }
//    }
//}
///**
//* @brief   function to write offset for each cell
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  num_cell   number of cells
//*/
//void VtkWriterManager::WirteGridCellOffsets(FILE* fp,
//    const bool bool_binary, const DefUint dims,
//    const DefLUint num_cell, OutputDataFormat& output_data_format) {
//
//    DefUint num_offsets = 4;
//    if (dims == 2) {
//        num_offsets = 4;  // VTK_PIXEL
//    } else if (dims == 3) {
//        num_offsets = 8;  // VTK_VOXEL
//    }
//    if (bool_binary) {
//        std::vector<uint8_t> vec_uint8{}, vec_base64{};
//        for (DefLUint i = 0; i < num_cell; ++i) {
//            base64_instance_.add_to_vec_char(
//                output_data_format.output_uint_.cast_type(
//                    i * num_offsets), &vec_uint8);
//        }
//        base64_instance_.Encode(&vec_uint8, &vec_base64);
//        for (const auto& iter : vec_base64) {
//            fprintf_s(fp, "%c", iter);
//        }
//        fprintf_s(fp, "\n");
//    } else {
//        std::string str_format = "   "
//            + output_data_format.output_uint_.printf_format_ + "\n";
//        for (DefLUint i = 0; i < num_cell; ++i) {
//            fprintf_s(fp, str_format.c_str(), i * num_offsets);
//        }
//    }
//}
///**
//* @brief   function to write cell type.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  num_cell   number of cells.
//*/
//void VtkWriterManager::WirteGridCellTypes(FILE* fp,
//    bool bool_binary, const DefUint dims, const DefLUint num_cell) {
//
//    std::uint8_t cell_type = 8;
//    if (dims == 2) {
//        cell_type = 8;  // VTK_PIXEL
//    } else if (dims == 3) {
//        cell_type = 11;  // VTK_VOXEL
//    }
//    if (bool_binary) {
//        std::vector<uint8_t> vec_uint8{}, vec_base64{};
//        for (DefLUint i = 0; i < num_cell; ++i) {
//            base64_instance_.add_to_vec_char(cell_type, &vec_uint8);
//        }
//        base64_instance_.Encode(&vec_uint8, &vec_base64);
//        for (const auto& iter : vec_base64) {
//            fprintf_s(fp, "%c", iter);
//        }
//        fprintf_s(fp, "\n");
//    } else {
//        for (DefLUint i = 0; i < num_cell; ++i) {
//            fprintf_s(fp, "   %u\n", cell_type);
//        }
//    }
//}
///**
//* @brief   function to write pieces of geometries in vtu.
//* @param[in]  fp   vertexer vertex to output file.
//*/
//void VtkWriterManager::WriteGeometry(FILE* fp,
//    bool bool_binary, OutputDataFormat& output_data_format,
//    std::shared_ptr<criterion::CriterionManager> ptr_criterion_manager) {
//
//    std::string str_temp;
//
//    for (const auto& iter : ptr_criterion_manager->vec_ptr_geometries_) {
//        DefSizet cell_number = CalculateNumOfGeometryCells(*iter);
//        str_temp.assign("  <Piece NumberOfvertexs=\""
//            + std::to_string(iter->vec_coordinate_origin_.size())
//            + "\" NumberOfCells=\"" + std::to_string(cell_number)
//            + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//        fprintf_s(fp, "   <vertexs>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_real_.format_name_
//            + "\" NumberOfComponents=\"3\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write vertex data
//        WirteGeometryCoordinates(fp, bool_binary,
//            ptr_criterion_manager->k0GeoDims_, output_data_format, *iter);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        fprintf_s(fp, "   </vertexs>\n");
//        fprintf_s(fp, "   <Cells>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_uint_.format_name_
//            + "\" Name=\"connectivity\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell connectiveity
//        WirteGeometryCellConnectivity(
//            fp, bool_binary, output_data_format, *iter);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        str_temp.assign("    <DataArray type=\""
//            + output_data_format.output_uint_.format_name_
//            + "\" Name=\"offsets\" format=\""
//            + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell offsets
//        WirteGeometryCellOffsets(
//            fp, bool_binary, output_data_format, *iter);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        str_temp.assign("    <DataArray type=\"UInt8\" Name=\"types\" "
//            "format=\"" + k0StrVtkAsciiOrBinary_ + "\">\n");
//        fprintf_s(fp, str_temp.c_str());
//
//        // write cell types
//        WirteGeometryCellTypes(
//            fp, bool_binary, *iter);
//
//        fprintf_s(fp, "    </DataArray>\n");
//        fprintf_s(fp, "   </Cells>\n");
//        fprintf_s(fp, "  </Piece>\n");
//    }
//}
///**
//* @brief   function to calculate number of cells.
//*          according to geometry_cell_type_.
//* @param[in]  geo_info   an instance of geometry information.
//* @return  number of cells.
//*/
//DefSizet VtkWriterManager::CalculateNumOfGeometryCells(
//    const criterion::GeometryInfoInterface& geo_info) {
//    DefSizet num_cell = 0;
//    switch (geo_info.geometry_cell_type_) {
//    case   criterion::EGeometryCellType::kPolyLine:
//        num_cell = geo_info.vec_coordinate_origin_.size() - 1;
//        break;
//    default:
//        int rank_id = 0;
//#ifdef ENABLE_MPI
//        rank_id = mpi::MpiManager::GetInstance()->rank_id_;
//#endif  // ENABLE_MPI
//        io::LogError("GeometryOutputType is undefined for calculing cell"
//            " numbers in CalculateGemetryCells.");
//        break;
//    }
//    return num_cell;
//}
///**
//* @brief   function to write geometry coordinates.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information.
//*/
//void VtkWriterManager::WirteGeometryCoordinates(FILE* fp,
//    const bool bool_binary, const DefUint dims,
//    OutputDataFormat& output_data_format,
//    const criterion::GeometryInfoInterface& geo_info) {
//
//    std::string str_format = "   "
//        + output_data_format.output_real_.printf_format_ + " "
//        + output_data_format.output_real_.printf_format_ + " "
//        + output_data_format.output_real_.printf_format_ + "\n";
//    DefReal x_offset = geo_info.k0RealOffest_.at(kXIndex),
//        y_offset = geo_info.k0RealOffest_.at(kYIndex);
//
//    if (dims == 2) {
//        if (bool_binary) {
//            std::vector<uint8_t> vec_uint8{}, vec_base64{};
//            for (const auto& iter : geo_info.vec_coordinate_origin_) {
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(
//                        iter.coordinate.at(kXIndex) - x_offset), &vec_uint8);
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(
//                        iter.coordinate.at(kYIndex) - y_offset), &vec_uint8);
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(0.), &vec_uint8);
//            }
//            base64_instance_.Encode(&vec_uint8, &vec_base64);
//            for (const auto& iter : vec_base64) {
//                fprintf_s(fp, "%c", iter);
//            }
//            fprintf_s(fp, "\n");
//        } else {
//            for (const auto& iter : geo_info.vec_coordinate_origin_) {
//                fprintf_s(fp, str_format.c_str(),
//                    iter.coordinate.at(kXIndex) - x_offset,
//                    iter.coordinate.at(kYIndex) - y_offset,
//                    0.);
//            }
//        }
//    } else if (dims == 3) {
//        DefReal z_offset = geo_info.k0RealOffest_.at(kZIndex);
//        if (bool_binary) {
//            std::vector<uint8_t> vec_uint8{}, vec_base64{};
//            for (const auto& iter : geo_info.vec_coordinate_origin_) {
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(
//                        iter.coordinate.at(kXIndex) - x_offset), &vec_uint8);
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(
//                        iter.coordinate.at(kYIndex) - y_offset), &vec_uint8);
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_real_.cast_type(
//                        iter.coordinate.at(kZIndex) - z_offset), &vec_uint8);
//            }
//            base64_instance_.Encode(&vec_uint8, &vec_base64);
//            for (const auto& iter : vec_base64) {
//                fprintf_s(fp, "%c", iter);
//            }
//            fprintf_s(fp, "\n");
//        } else {
//            for (const auto& iter : geo_info.vec_coordinate_origin_) {
//                fprintf_s(fp, str_format.c_str(),
//                    iter.coordinate.at(kXIndex) - x_offset,
//                    iter.coordinate.at(kYIndex) - y_offset,
//                    iter.coordinate.at(kZIndex) - z_offset);
//            }
//        }
//    }
//}
///**
//* @brief   function to write connectivity relation of cells.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information.
//*/
//void VtkWriterManager::WirteGeometryCellConnectivity(FILE* fp,
//    const bool bool_binary, OutputDataFormat& output_data_format,
//    const criterion::GeometryInfoInterface& geo_info) {
//
//    switch (geo_info.geometry_cell_type_) {
//    case   criterion::EGeometryCellType::kPolyLine: {
//        DefLUint i_vertex = 0;
//        if (bool_binary) {
//            std::vector<uint8_t> vec_uint8{}, vec_base64{};
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_uint_.cast_type(
//                        i_vertex), &vec_uint8);
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_uint_.cast_type(
//                        i_vertex + 1), &vec_uint8);
//                ++i_vertex;
//            }
//            base64_instance_.Encode(&vec_uint8, &vec_base64);
//            for (const auto& iter : vec_base64) {
//                fprintf_s(fp, "%c", iter);
//            }
//            fprintf_s(fp, "\n");
//        } else {
//            std::string str_format = "   "
//                + output_data_format.output_uint_.printf_format_ + " "
//                + output_data_format.output_uint_.printf_format_ + "\n";
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                fprintf_s(fp, str_format.c_str(), i_vertex,
//                    i_vertex + 1);
//                ++i_vertex;
//            }
//        }
//        break;
//    }
//    default:
//        break;
//    }
//}
///**
//* @brief   function to write offset for each cell.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information.
//*/
//void VtkWriterManager::WirteGeometryCellOffsets(FILE* fp,
//    const bool bool_binary, OutputDataFormat& output_data_format,
//    const criterion::GeometryInfoInterface& geo_info) {
//
//    switch (geo_info.geometry_cell_type_) {
//    case   criterion::EGeometryCellType::kPolyLine: {
//        DefLUint i_vertex = 1;
//        if (bool_binary) {
//            std::vector<uint8_t> vec_uint8{}, vec_base64{};
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                base64_instance_.add_to_vec_char(
//                    output_data_format.output_uint_.cast_type(
//                        i_vertex * 2), &vec_uint8);
//                ++i_vertex;
//            }
//            base64_instance_.Encode(&vec_uint8, &vec_base64);
//            for (const auto& iter : vec_base64) {
//                fprintf_s(fp, "%c", iter);
//            }
//            fprintf_s(fp, "\n");
//        } else {
//            std::string str_format = "   "
//                + output_data_format.output_uint_.printf_format_ + "\n";
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                fprintf_s(fp, str_format.c_str(), i_vertex * 2);
//                ++i_vertex;
//            }
//        }
//        break;
//    }
//    default:
//        break;
//    }
//}
///**
//* @brief   function to write cell type.
//* @param[in]  fp   vertexer vertex to output file.
//* @param[in]  geo_info   an instance of geometry information.
//*/
//void VtkWriterManager::WirteGeometryCellTypes(FILE* fp,
//    const bool bool_binary, const criterion::GeometryInfoInterface& geo_info) {
//
//    switch (geo_info.geometry_cell_type_) {
//    case   criterion::EGeometryCellType::kPolyLine:
//        if (bool_binary) {
//            std::vector<uint8_t> vec_uint8{}, vec_base64{};
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                base64_instance_.add_to_vec_char(static_cast<std::uint8_t>(
//                    geo_info.geometry_cell_type_), &vec_uint8);
//            }
//            base64_instance_.Encode(&vec_uint8, &vec_base64);
//            for (const auto& iter : vec_base64) {
//                fprintf_s(fp, "%c", iter);
//            }
//            fprintf_s(fp, "\n");
//        } else {
//            for (auto iter = geo_info.vec_coordinate_origin_.begin();
//                iter != geo_info.vec_coordinate_origin_.end() - 1; ++iter) {
//                fprintf_s(fp, "   %u\n",
//                    static_cast<std::uint8_t>(geo_info.geometry_cell_type_));
//            }
//        }
//        break;
//    default:
//        break;
//    }
//}
//}  // end namespace io
//}  // end amrproject
//}  // end namespace rootproject
