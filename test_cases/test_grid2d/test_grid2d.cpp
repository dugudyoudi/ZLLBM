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
class GridInfoTest :public GridInfoInterface {
 public:
    void SetNumberOfVecElements() override {};
    void InitialGridNode(const DefSFBitset& bitset_in) override {};
    GridInfoTest() {
        this->node_type_ = "GridInfoTest";
    }
};
class GridInfoTestCreator :
    public GridInfoCreatorInterface {
 public:
    std::shared_ptr<GridInfoInterface>
        CreateGridInfo()override {
        return std::make_shared<GridInfoTest>();
    };
};
class TrackingGridInfoTest :public TrackingGridInfoInterface {
 public:
    TrackingGridInfoTest() {
        this->node_type_ = "TrackingGridInfoTest";
    }
};
class TrackingGridInfoTestCreator :
    public TrackingGridInfoCreatorInterface {
 public:
    std::shared_ptr<TrackingGridInfoInterface>
        CreateTrackingGridInfo()override {
        return std::make_shared<TrackingGridInfoTest>();
    };
};
class SolverTest :public SolverInterface {
 public:
    std::string GetSolverMethod() override {
        return "A temporal solver for test.";
    };
    void SolverInitial() override {};
    void RunSolver() override {};
    SolverTest() {
        ptr_grid_info_creator_ = std::make_unique<GridInfoTestCreator>();
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
    DefSizet max_refinement_level = 2;  // maximum refinement level

    // instantiate modules will be used
    AmrManager* amr_instance = AmrManager::GetInstance();

    // must call DefaultInitialize first for mpi initialization
    amr_instance->DefaultInitialization(
        dims, max_refinement_level, argc, argv);

    amr_instance->ptr_io_manager_->bool_binary_ = false;
    amr_instance->ptr_io_manager_->vtk_instance_.vtk_ghost_cell_option_ =
        amrproject::EVtkWriterGhostCellOption::kPartitionMultiBlock;
    //amr_instance->ptr_io_manager_->bool_binary_ = true;

    // grid related parameters //
    GridManager2D* grid_manager = dynamic_cast<GridManager2D*>(
        amr_instance->ptr_grid_manager_.get());
    grid_manager->vec_ptr_tracking_info_creator.push_back(
        std::make_unique<TrackingGridInfoTestCreator>());
    // domain size
    grid_manager->k0DomainSize_.at(kXIndex) = 2.;
    grid_manager->k0DomainSize_.at(kYIndex) = 2.;
    // grid space
    grid_manager->k0DomainDx_.at(kXIndex) = 0.02;
    // end grid related parameters //

    // geometry related parameters //
    GeometryInfoConnection2DCreator geo_creator;
    amr_instance->ptr_criterion_manager_->vec_ptr_geometries_.push_back(
        geo_creator.CreateGeometryInfo());
    // convert base ptr to derived ptr for visiting members in derived class
    GeometryInfoConnection2D* ptr_geo_temp =
        dynamic_cast<GeometryInfoConnection2D*>(amr_instance->
        ptr_criterion_manager_->vec_ptr_geometries_.at(0).get());
    ptr_geo_temp->ptr_tracking_grid_info_creator_ = grid_manager->vec_ptr_tracking_info_creator.at(0).get();
    ptr_geo_temp->geometry_center_ = { 1, 1 };
    ptr_geo_temp->k0XIntExtendPositive_ = {3, 6 };
    ptr_geo_temp->k0XIntExtendNegative_ =
        std::vector<DefUint>(max_refinement_level + 1 , 2);
    ptr_geo_temp->k0YIntExtendPositive_ = { 2, 4 };
    ptr_geo_temp->k0YIntExtendNegative_ = { 2, 4 };
    ptr_geo_temp->k0IntInnerExtend_ = { 2, 2 };
    ///< number of extened layers

    ptr_geo_temp->grid_extend_type_ =
        EGridExtendType::kInAndOut;
    ptr_geo_temp->bool_vec_forces_ = false;
    ptr_geo_temp->i_level_ = max_refinement_level;
    /* used for generating predefined geometries, number of input parameters
    is based on the type of geometry_shape_*/
    DefReal dx = grid_manager->k0DomainDx_.at(kXIndex)
        / std::pow(2, max_refinement_level);
    ptr_geo_temp->InitialGeometry(dx, DefaultGeoShapeType::kCircle,
        *(amr_instance->ptr_criterion_manager_->ptr_default_geo_manager_));
    std::vector<DefReal> coordi_min, coordi_max;
    ptr_geo_temp->InitialConnection(&coordi_min, &coordi_max);
    // end geometry related parameters //

    amr_instance->SetupParameters();

    // set grid node type and solver at all levels the same
    std::shared_ptr<SolverCreatorTest> ptr_solver_creator =
        std::make_shared<SolverCreatorTest>();
    amr_instance->SetTheSameLevelDependentInfoForAllLevels(
        ptr_solver_creator.get());

    amr_instance->InitializeSimulation();

#ifdef ENABLE_MPI
    int myid, numprocs, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Get_processor_name(processor_name, &namelen);

    // if (myid == 0) printf("number of processes: %d\n", numprocs);
    printf("Mpi test (test_grid.cpp) finished on %s: node %d \n",
        processor_name, myid);
#else
    printf("Serial test (test_grid2d.cpp) finished");
#endif

    amr_instance->FinalizeSimulation();

}
