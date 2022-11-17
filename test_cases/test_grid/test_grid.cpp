//    Copyright (c) 2022, Zhengliang Liu
//    All rights reserved

/**
* @file test_grid.cpp
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  test grid generation and related MPI functions
*/
#include "amr_manager.h"
#include "criterion/geometry_info_connection.h"
using namespace rootproject;
using namespace rootproject::amrproject;
class GridInfoTest :public grid::GridInfoInterface {
public:
    void set_number_of_vec_elements() override {};
    void InitialGridNode(const DefSFBitset& bitset_in) override {};
    GridInfoTest() {
        this->node_type_ = "GridInfoTest";
    }
};
class GridInfoTestCreator :
    public grid::GridInfoCreatorInterface {
public:
    std::shared_ptr<grid::GridInfoInterface>
        CreateGridInfo()override {
        return std::make_shared<GridInfoTest>();
    };
};
class SolverTest :public SolverInterface {
 public:
    std::string GetSolverMethod() override {
        return "A temperal sovler for test.";
    };
    void SolverInitial() override {};
    void RunSolver() override {};
    SolverTest() {
        ptr_grid_info_creator_ = std::make_shared<GridInfoTestCreator>();
    }
};
class SolverCreatorTest :public SolverCreatorInterface {
    std::shared_ptr<SolverInterface> CreateSolver() override {
        std::shared_ptr<SolverTest> ptr_temp =
            std::make_shared<SolverTest>();
        return ptr_temp;
    }
};
int main(int argc, char* argv[]) {
    DefUint dims = 2;   // dimension
    DefSizet max_refinement_level = 1;  // maximum refinement level


    std::shared_ptr<criterion::GeometryInfoConnection2D>
        ptr_geo_info_ = std::make_shared<criterion::GeometryInfoConnection2D>();
    DefReal ds_origin = 0.1;
    std::array<DefReal, 2> center = { 4, 4 };
    DefSizet num_point;
    DefReal r = 0.5;

    ptr_geo_info_->geometry_cell_type_
        = criterion::EGeometryCellType::kTriangle;
    ptr_geo_info_->SetIndex();
    ptr_geo_info_->SetupConnectionParameters(
        criterion::EGeometryCellType::kTriangle);
    ptr_geo_info_->geometry_center_ =
    { center[kXIndex],  center[kYIndex] + r * 2 };
    num_point = 5;
    ptr_geo_info_->coordinate_origin_ =
        std::vector<criterion::GeometryCoordinate2D>(
            num_point, { 0., 0. });
    std::vector<DefReal> xcor{ -0.2, 0.3, 0.2, -0.2, 0 };
    std::vector<DefReal> ycor{ 0, 0, 0.4, 0.4, 0.3 };
    for (DefUint i = 0; i < num_point; ++i) {
        ptr_geo_info_->coordinate_origin_.at(i).coordinate.at(kXIndex) =
            xcor[i] + ptr_geo_info_->geometry_center_[kXIndex];
        ptr_geo_info_->coordinate_origin_.at(i).coordinate.at(kYIndex) =
            ycor[i] * ds_origin + ptr_geo_info_->geometry_center_[kYIndex];
    }
    ptr_geo_info_->connection_relation_.push_back({ 0,1,4 });
    ptr_geo_info_->connection_relation_.push_back({ 1,2,4 });
    ptr_geo_info_->connection_relation_.push_back({ 2,3,4 });
    ptr_geo_info_->connection_relation_.push_back({ 3,0,4 });
    DefReal ds_max = ds_origin / 3.;
    DefSizet i_input_level = 0;
    DefSizet num_vertex_origin = ptr_geo_info_
        ->coordinate_origin_.size();
    ptr_geo_info_->InitialConnection();
    std::set<std::pair<std::pair<DefSizet, DefSizet>, std::pair<DefSizet,
        DefSizet>>> edges_for_bisect, edges_remain;
    edges_for_bisect.insert({ {0, 4}, {0, 0} });
    edges_for_bisect.insert({ {0, 1}, {0, 0} });
    edges_for_bisect.insert({ {0, 4}, {0, 1} });
    edges_for_bisect.insert({ {0, 2}, {0, 1} });
    edges_for_bisect.insert({ {0, 4}, {0, 2} });
    ptr_geo_info_->BisectEdgeOnce(i_input_level, ds_max,
        edges_for_bisect, &edges_remain);

    DefReal ds_min = ds_origin * 3.;
    edges_remain.clear();
    std::set<std::pair<std::pair<DefSizet, DefSizet>, std::pair<DefSizet,
        DefSizet>>> edges_for_merge;
    std::pair<DefSizet, DefSizet> vertex0, vertex1;
    DefSizet num_vertex = ptr_geo_info_->
        connection_vertex_given_level_.size();
    DefSizet num_edge = ptr_geo_info_->
        connection_edge_given_level_.size();
    DefSizet num_surface = ptr_geo_info_->
        connection_surface_given_level_.size();
    std::array<DefReal, 2> coorid_mid;
    coorid_mid[kXIndex] = (ptr_geo_info_->coordinate_origin_.at(0)
        .coordinate[kXIndex] + ptr_geo_info_->coordinate_origin_
        .at(1).coordinate[kXIndex]) / 2;
    coorid_mid[kYIndex] = (ptr_geo_info_->coordinate_origin_.at(0)
        .coordinate[kYIndex] + ptr_geo_info_->coordinate_origin_
        .at(1).coordinate[kYIndex]) / 2;


    vertex0 = { 1, 0 }; vertex1 = { 0, 4 };
    edges_for_merge.insert({ vertex0 , vertex1 });
    ptr_geo_info_->MergeEdgeOnce(i_input_level, ds_min,
        edges_for_merge, &edges_remain);

    DefSizet  i_point = 0;
    edges_for_merge.clear();
    for (const auto& iter : ptr_geo_info_->vertex_given_level_
        .at(i_input_level + 1).vec_vertex_cooridinate) {
        vertex0 = { i_input_level + 1,  i_point };
        for (const auto& iter_vertex : iter.map_linked_vertices_level
            .at(i_input_level)) {
            if (vertex0 > iter_vertex) {
                edges_for_merge.insert({ vertex0, iter_vertex });
            } else {
                edges_for_merge.insert({ iter_vertex, vertex0 });
            }
            //std::cout << vertex0.first << " " << vertex0.second << "; " <<
            //    iter_vertex.first << " " << iter_vertex.second << std::endl;
        }
        ++i_point;
    }

    ptr_geo_info_->MergeEdgeOnce(i_input_level, ds_min,
        edges_for_merge, &edges_remain);


    // instantiate modules will be used
//    AmrManager* amr_instance = AmrManager::GetInstance();
//    DefInt mpi_id = 0;
//#ifdef ENABLE_MPI
//    using rootproject::mpi::MpiManager;
//    MpiManager* mpi_instance = MpiManager::GetInstance();
//#endif
//    // must call DefaultInitialize fisrt for mpi initialization
//    amr_instance->DefaultInitialization(
//        dims, max_refinement_level, argc, argv);
//
//    // grid related parameters //
//    // domain size
//    amr_instance->ptr_grid_manager_->k0DomainSize_.at(kXIndex) = 2.;
//    amr_instance->ptr_grid_manager_->k0DomainSize_.at(kYIndex) = 2.;
//    // grid space
//    amr_instance->ptr_grid_manager_->k0DomainDx_.at(kXIndex) = 0.02;
//    // end grid related parameters //

    //// geometry related parameters //
    //amr_instance->ptr_criterion_manager_->vec_ptr_geometries_.push_back(
    //std::make_shared<criterion::GeometryInfoStl>());
    //// convert base ptr to derived ptr for visiting members in derived class
    //std::shared_ptr<criterion::GeometryInfoStl> geo_ptr_temp =
    //    std::dynamic_pointer_cast<criterion::GeometryInfoStl>
    //    (amr_instance->ptr_criterion_manager_->vec_ptr_geometries_.at(0));
    //geo_ptr_temp->geometry_center_ = { 1, 1 };
    //geo_ptr_temp->k0DefaultGeoShapeType_ = criterion::DefaultGeoShapeType::kCircle;
    ////geo_ptr_temp->k0GeoGhostType_ =
    ////    geometry::GeometryGhostType::kResetExtendLayers;
    //geo_ptr_temp->bool_vec_forces_ = false;
    ///* used for generating predefined geometries, number of input paramters
    //is based on the type of geometry_shape_*/
    //geo_ptr_temp->geometry_input_parameters_.push_back(0.5);
    //// end geometry related parameters //

    //// set grid node type and sovler at all levels the same
    //std::shared_ptr<SolverCreatorTest> ptr_sovler_creator =
    //    std::make_shared<SolverCreatorTest>();
    //amr_instance->SetTheSameLevelDependentInfoForAllLevels(
    //    ptr_sovler_creator);

    //amr_instance->InitialSimulation();

    //amr_instance->FinalizeSimulation();
#ifdef ENABLE_MPI
    int myid, numprocs, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    // MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Get_processor_name(processor_name, &namelen);

    // if (myid == 0) printf("number of processes: %d\n", numprocs);
    printf("Mpi test (test_grid.cpp) finished on %s: node %d \n",
        processor_name, myid);

    MPI_Finalize();
#else
    printf("Serial test (test_grid.cpp) finished");
#endif
}
