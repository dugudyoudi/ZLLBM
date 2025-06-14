//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file auxiliary_inline_func.h
* @author Zhengliang Liu
* @brief amp project inline functions
* @date  2022-9-3
* @note .
*/
#ifndef SOURCE_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#define SOURCE_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#include <array>
#include <limits>
#include "../defs_libs.h"
namespace rootproject {
inline DefReal Square(DefReal x) {
    return x * x;
}
/**
* @brief function to calculate 2 power N.
* @param[in]  n        power exponent.
*/
inline DefSizet TwoPowerN(DefSizet n) {
    DefSizet two = 1;
    two <<= n;
    return two;
}
inline DefAmrLUint TwoPowerN(DefInt n) {
    DefAmrLUint two = 1;
    two <<= n;
    return two;
}
constexpr DefReal SqrtNewtonRaphson(DefReal x, DefReal current, DefReal previous) {
    return current == previous ? current : SqrtNewtonRaphson(x, 0.5 * (current + x / current), current);
}
constexpr DefReal SqrtConstexpr(DefReal x) {
    return x >= 0 && x < std::numeric_limits<DefReal>::infinity()
        ? SqrtNewtonRaphson(x, x, 0)
        : std::numeric_limits<DefReal>::quiet_NaN();
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
