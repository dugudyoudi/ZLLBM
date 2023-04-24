//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file auxiliary_inline_func.h
* @author Zhengliang Liu
* @brief amp project inline functions
* @date  2022-9-3
* @note .
*/
#ifndef ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#define ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#include <array>
#include <fstream>
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
inline DefReal Square(DefReal x) {
    return x * x;
}
/**
* @brief function to calculate 2 power N.
* @param[in]  n        power exponent.
*/
inline uint64_t TwoPowerN(uint64_t n) {
    uint64_t two = 1;
    two <<= n;
    return two;
}
inline uint32_t TwoPowerN(uint32_t n) {
    uint32_t two = 1;
    two <<= n;
    return two;
}
inline uint16_t TwoPowerN(uint16_t n) {
    uint16_t two = 1;
    two <<= n;
    return two;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_