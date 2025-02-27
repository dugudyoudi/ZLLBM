//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.h
* @author Zhengliang Liu
* @brief define classes to manage immersed boundary method.
* @date  2024-2-03
*/
#ifndef SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
#define SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
#include <unordered_map>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include "criterion/geometry_info_origin.h"
#include "./lbm_interface.h"
#include "immersed_boundary/geometry_ib_shape.h"
namespace rootproject {
namespace lbmproject {
class GeoIBShapeReader;
class GeometryVertexImmersedBoundary : public amrproject::GeometryVertex {
 public:
    std::array<DefReal, 3> velocity_{};
    std::array<DefReal, 3> force_{};
    DefReal area_;
};
class FsiImmersedBoundary {
 public:
    void CalculateBodyForce(const GridInfoLbmInteface& grid_info,
        std::unordered_map<DefSizet, std::unique_ptr<amrproject::GeometryVertex>>* const ptr_geo_vertices,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes);
    void ClearNodesRecordForIB();
    void SetStencilDis(const DefInt stencil_dis) { stencil_dis_ = stencil_dis;}
    void SetWriteIBForce(const bool write_ib_force) { write_ib_force_ = write_ib_force;}

 protected:
    bool write_ib_force_ = false;
    DefInt stencil_dis_ = 2;
    DefReal ib_stiffness_ = 2.;
    DefInt k0IndexForceX_ = 0, k0IndexForceY_ = k0IndexForceX_ + 1, k0IndexForceZ_ = k0IndexForceY_ + 1;

    void DirectForcingScheme2D(const DefReal dt_lbm, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const amrproject::DomainInfo& domain_info, const std::function<void(const DefReal, GridNodeLbm* const,
        DefReal* const, std::vector<DefReal>* const)>& func_compute_macro,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
        GeometryVertexImmersedBoundary* const ptr_geo_vertex) const;
    void DirectForcingScheme3D(const DefReal dt_lbm, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const amrproject::DomainInfo& domain_info, const std::function<void(const DefReal, GridNodeLbm* const,
        DefReal* const, std::vector<DefReal>* const)>& func_compute_macro,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
        GeometryVertexImmersedBoundary* const ptr_geo_vertex) const;
    void CopyIBNodeToBuffer2D(const GridNodeLbm& node, char* const ptr_node_buffer) const;
    void CopyIBNodeToBuffer3D(const GridNodeLbm& node, char* const ptr_node_buffer) const;
    void ReadIBNodeFromBuffer2D(const char* ptr_node_buffer, GridNodeLbm* const ptr_node) const;
    void ReadIBNodeFromBuffer3D(const char* ptr_node_buffer, GridNodeLbm* const ptr_node) const;

 private:
    static DefMap<DefInt> map_ib_node_for_reset_;

    DefReal StencilDisOne(DefReal dis) const;
    DefReal StencilDisTwo(DefReal dis) const;

    // mpi related
 public:
    void SendNReceiveNodesForImmersedBoundary(const DefInt dim, const DefInt i_level,
        const amrproject::MpiManager& mpi_manager,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes);

 protected:
    std::map<int, DefMap<DefInt>> map_ib_node_for_mpi_send_;
};

class GeometryInfoImmersedBoundary : public amrproject::GeometryInfoOrigin, public FsiImmersedBoundary {
 public:
    std::unique_ptr<amrproject::GeometryVertex> GeoInfoVertexCreator() const override {
        return std::make_unique<GeometryVertexImmersedBoundary>();
    }
    void SetupGeometryInfo(const DefReal time, const amrproject::MpiManager& mpi_manager,
        const amrproject::GridInfoInterface& sfbitset_aux) override;
    void ReadAndSetGeoParameters(const DefInt level,
        std::map<std::string, amrproject::ParserData>* const ptr_geo_parameters) override;
    void WriteTimeHisLagrangianForce(const DefReal time, const DefReal dx_background) const;
    explicit GeometryInfoImmersedBoundary(const DefInt dims) : amrproject::GeometryInfoOrigin(dims) {
        this->node_type_ = "GeometryInfoImmersedBoundary";
    }
};
class GeometryInfoImmersedBoundaryCreator :public amrproject::GeometryInfoOriginCreator {
 public:
    std::shared_ptr<amrproject::GeometryInfoInterface> CreateGeometryInfo(const DefInt dims) const override {
        std::shared_ptr<GeometryInfoImmersedBoundary> ptr_tmp = std::make_shared<GeometryInfoImmersedBoundary>(dims);
        return ptr_tmp;
    };
};
class GeoIBTypeReader : public amrproject::GeoTypeReader {
 public:
    std::shared_ptr<amrproject::GeometryInfoInterface> ReadGeoType(
        const DefInt dims, const std::string& geo_type = "origin") const override;
    std::unique_ptr<amrproject::GeoShapeReader> CreateShapeReader() const override {
        return std::make_unique<GeoIBShapeReader>();
    }
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
