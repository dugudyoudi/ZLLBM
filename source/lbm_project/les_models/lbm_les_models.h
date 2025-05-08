//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_les_models.h
* @author Zhengliang Liu
* @brief define classes to manage models for large eddy simulation.
* @date  2025-2-20
*/
#ifndef SOURCE_LBM_PROJECT_LES_MODELS_LBM_LES_MODELS_H_
#define SOURCE_LBM_PROJECT_LES_MODELS_LBM_LES_MODELS_H_
#include <vector>
#include "../defs_libs.h"
namespace rootproject {
namespace lbmproject {
struct GridNodeLbm;
class SolverLbmInterface;
/**
* @brief enumerate LES model type
*/
enum class ELbmLesModelType {
    kSmagorinsky,
    kUndefined
};
/**
* @brief interface class to manage LES model
*/
class LesModelInterface {
 public:
    virtual DefReal CalSgsRelaxationTimeWithoutForce(const DefReal dt_lbm,
        const DefReal tau_0, const std::vector<DefReal>& feq,
        const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const = 0;
    virtual DefReal CalSgsRelaxationTimeWithForce(const DefReal dt_lbm,
        const DefReal tau_0, const std::vector<DefReal>& feq, const std::vector<DefReal>& force,
        const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const = 0;
    virtual ~LesModelInterface() {}
};
/**
* @brief interface class to manage LES model
*/
class LesModelSmagorinsky : public LesModelInterface{
    DefReal CalSgsRelaxationTimeWithoutForce(const DefReal dt_lbm,
        const DefReal tau_0, const std::vector<DefReal>& feq,
        const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const override;
    DefReal CalSgsRelaxationTimeWithForce(const DefReal dt_lbm,
        const DefReal tau_0, const std::vector<DefReal>& feq,  const std::vector<DefReal>& force,
        const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const override;

 protected:
    DefReal k0Smagorinsky = 18.*0.16*0.16;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_LES_MODELS_LBM_LES_MODELS_H_

