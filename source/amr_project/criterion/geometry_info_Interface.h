//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define classes to store geometry information.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#define SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
#include <vector>
#include <array>
#include <memory>
#include <unordered_map>
#include <map>
#include <string>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
#include "criterion/geometry_default_shape.h"
namespace rootproject {
namespace amrproject {
class TrackingGridInfoCreatorInterface;
class GhostGridInfoCreatorInterface;
class SFBitsetAuxInterface;
class GridInfoInterface;
class GeoShapeInterface;
class GeoShapeReader;
class MpiManager;
/**
* @struct GeometryVertex
* @brief struct used to store information of a geometry vertex
*/
struct GeometryVertex {
    std::array<DefReal, 3> coordinate_{};
    virtual ~GeometryVertex() {}
};
/**
* @class GeometryInfoInterface
* @brief class used to store information of a geometry
*/
class GeometryInfoInterface {
 protected:
    // information of geometry itself
    DefInt k0GeoDim_ = 0;
    DefInt computational_cost_ = 1;
    DefInt i_level_ = 0;
    DefInt i_geo_ = ~0;
    EGeometryCellType geometry_cell_type_ = EGeometryCellType::kUndefined;
    EGridExtendType grid_extend_type_ = EGridExtendType::kSameInAllDirections;
    EGeometryStatus geometry_status_ = EGeometryStatus::kVirtual;
    std::string name_, node_type_;
    bool need_update_shape_ = false;

    TrackingGridInfoCreatorInterface* ptr_tracking_grid_info_creator_ = nullptr;

    // number of extended layer based on geometry
    std::vector<DefInt> k0XIntExtendPositive_, k0XIntExtendNegative_, k0YIntExtendPositive_,
        k0YIntExtendNegative_, k0ZIntExtendPositive_, k0ZIntExtendNegative_;
        ///< number of extened layers

    /* number of layer extended inside the geometry
     at (i_level - 1) refinement level*/
    std::vector<DefInt> k0IntInnerExtend_;
    ///< number of extened layers inside the geometry

    std::array<DefReal, 3> geometry_center_{};
    std::array<DefReal, 3> flood_fill_origin_{};
    std::array<DefReal, 3> k0RealMin_{};

 public:
    // get and set functions
    bool GetNeedUpdateShape() const {
        return need_update_shape_;
    }
    DefInt GetComputationalCost() const {
        return computational_cost_;
    }
    DefInt GetGeoIndex() const {
        return i_geo_;
    }
    DefInt GetLevel() const {
        return i_level_;
    }
    EGeometryCellType GetCellType() const {
        return geometry_cell_type_;
    }
    EGridExtendType GetGridExtendType() const {
        return grid_extend_type_;
    }
    EGeometryStatus GetStatus() const {
        return geometry_status_;
    }
    std::string GetName() const {
        return name_;
    }
    TrackingGridInfoCreatorInterface* GetPtrTrackingGridInfoCreatorInterface() const {
        return ptr_tracking_grid_info_creator_;
    }
    std::vector<DefInt> GetXExtendPositive() const {
        return k0XIntExtendPositive_;
    }
    std::vector<DefInt> GetXExtendNegative() const {
        return k0XIntExtendNegative_;
    }
    std::vector<DefInt> GetYExtendPositive() const {
        return k0YIntExtendPositive_;
    }
    std::vector<DefInt> GetYExtendNegative() const {
        return k0YIntExtendNegative_;
    }
    std::vector<DefInt> GetZExtendPositive() const {
        return k0ZIntExtendPositive_;
    }
    std::vector<DefInt> GetZExtendNegative() const {
        return k0ZIntExtendNegative_;
    }
    std::vector<DefInt> GetInnerExtend() const {
        return k0IntInnerExtend_;
    }
    std::array<DefReal, 3> GetGeometryCenter() const {
        return geometry_center_;
    }
    std::array<DefReal, 3> GetFloodFillOrigin() const {
        return flood_fill_origin_;
    }
    std::array<DefReal, 3> GetOffset() const {
        return k0RealMin_;
    }
    virtual DefSizet GetNumOfGeometryPoints() const {
        return vec_vertices_.size();
    }

    void SetNeedUpdateShape(const bool need_update) {
        need_update_shape_ = need_update;
    }
    void SetComputationalCost(const DefInt computational_cost) {
        computational_cost_ = computational_cost;
    }
    void SetGeoIndex(const DefInt i_geo) {
        i_geo_ = i_geo;
    }
    void SetLevel(const DefInt i_level) {
        i_level_ = i_level;
    }
    void SetCellType(const EGeometryCellType geometry_cell_type) {
        geometry_cell_type_ = geometry_cell_type;
    }
    void SetGridExtendType(const EGridExtendType grid_extend_type) {
        grid_extend_type_ = grid_extend_type;
    }
    void SetStatus(const EGeometryStatus geometry_status) {
        geometry_status_ = geometry_status;
    }
    void SetName(const std::string& geo_name) {
        name_ = geo_name;
    }
    void SetPtrTrackingGridInfoCreator(
        TrackingGridInfoCreatorInterface* const ptr_tracking_grid_info_creator) {
        ptr_tracking_grid_info_creator_ = ptr_tracking_grid_info_creator;
    }
    void SetXExtendPositive(const std::vector<DefInt>& vec_extend) {
        k0XIntExtendPositive_ = vec_extend;
    }
    void SetXExtendNegative(const std::vector<DefInt>& vec_extend) {
        k0XIntExtendNegative_ = vec_extend;
    }
    void SetYExtendPositive(const std::vector<DefInt>& vec_extend) {
        k0YIntExtendPositive_ = vec_extend;
    }
    void SetYExtendNegative(const std::vector<DefInt>& vec_extend) {
        k0YIntExtendNegative_ = vec_extend;
    }
    void SetZExtendPositive(const std::vector<DefInt>& vec_extend) {
        k0ZIntExtendPositive_ = vec_extend;
    }
    void SetZExtendNegative(const std::vector<DefInt>& vec_extend) {
        k0ZIntExtendNegative_ = vec_extend;
    }
    void SetInnerExtend(const std::vector<DefInt>& vec_extend) {
        k0IntInnerExtend_ = vec_extend;
    }
    void SetGeometryCenter(const std::array<DefReal, 3>& array_center) {
        geometry_center_ = array_center;
    }
    void SetFloodFillOrigin(const std::array<DefReal, 3>& array_origin) {
        flood_fill_origin_ = array_origin;
    }
    void SetOffset(const std::array<DefReal, 3>& array_offset) {
        k0RealMin_ = array_offset;
    }

    void ChooseGridExtendType(const std::string type_string);
    virtual void ReadAndSetGeoParameters(const DefInt level,
        const std::map<std::string, std::string>& geo_parameters);

    std::vector<std::unique_ptr<GeometryVertex>> vec_vertices_{};
    // only vertex information on current rank will be stored,
    // vertex coordinates on other ranks can be found in vec_vertices_
    std::unordered_map<DefSizet, std::unique_ptr<GeometryVertex>> map_vertices_info_{};
    std::unique_ptr<GeoShapeInterface> ptr_geo_shape_;
    virtual void InitialGeometry(const DefReal dx);
    virtual void UpdateGeometry(const DefReal sum_t);
    virtual void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) = 0;

    std::vector<DefReal> GetFloodFillOriginArrAsVec() const {
        if (k0GeoDim_ == 2) {
            return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex]};
        } else {
            return {flood_fill_origin_[kXIndex], flood_fill_origin_[kYIndex], flood_fill_origin_[kZIndex]};
        }
    }
    virtual std::unique_ptr<GeometryVertex> GeoIndexVertexCreator() const {
        return std::make_unique<GeometryVertex>();
    }
    virtual std::unique_ptr<GeometryVertex> GeoInfoVertexCreator() const {
        return std::make_unique<GeometryVertex>();
    }

    void InstantiateGeometryVertexInfo(
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max, const SFBitsetAuxInterface& sfbitset_aux);
    virtual void SetupGeometryInfo(const DefReal time,
        const MpiManager& mpi_manager, const GridInfoInterface& grid_info);

        explicit GeometryInfoInterface(const DefInt dims) : k0GeoDim_(dims) {
        this->node_type_ = "GeometryInfoInterface";
        this->geometry_cell_type_ = EGeometryCellType::kPolyLine;
    }
    virtual ~GeometryInfoInterface() {}

 private:
    GeometryInfoInterface();
};
/**
* @class GeometryInfoCreatorInterface
* @brief interface class used to create geometry instance
*/
class GeometryInfoCreatorInterface {
 public:
    virtual std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo(const DefInt dims) const = 0;
    virtual ~GeometryInfoCreatorInterface() {}
};
class GeoTypeReader {
 public:
    virtual std::shared_ptr<GeometryInfoInterface> ReadGeoType(
        const DefInt dims, const std::string& geo_type = "origin") const;
    virtual std::unique_ptr<GeoShapeReader> CreateShapeReader() const {return std::make_unique<GeoShapeReader>();}
    virtual ~GeoTypeReader() {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_INTERFACE_H_
