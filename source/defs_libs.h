//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file def_libs.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief contain definations and libraries will be used in all modules.
*/

#ifndef ROOTPROJECT_SOURCE_DEF_LIB_H_
#define ROOTPROJECT_SOURCE_DEF_LIB_H_
// libraries
#include <bitset>
#include <cstdint>
#include <unordered_map>
#include <iostream>
#include "./config.h"  // configuration file generated
namespace rootproject {
using DefReal = double;
using DefInt = int;
using DefLInt = int;
using DefUint = unsigned int;
using DefLUint = unsigned int;
using DefTypeUint = unsigned int;
using DefSFCodeToUint = uint64_t;
using DefSizet = size_t;

// data type for mpi communication
#define MPI_REAL_DATA_TYPE MPI_DOUBLE
#define MPI_INT_DATA_TYPE MPI_INT
#define MPI_LINT_DATA_TYPE MPI_INT
#define MPI_UINT_DATA_TYPE MPI_UNSIGNED
#define MPI_LUINT_DATA_TYPE MPI_UNSIGNED
#define MPI_TYPEUINT_DATA_TYPE MPI_UNSIGNED
#define MPI_CODEUINT_DATA_TYPE MPI_UINT64_T
#if SIZE_MAX == UCHAR_MAX
#define MPI_SIZET_DATA_TYPE MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MPI_SIZET_DATA_TYPE MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MPI_SIZET_DATA_TYPE MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MPI_SIZET_DATA_TYPE MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MPI_SIZET_DATA_TYPE MPI_UNSIGNED_LONG_LONG
#else
#error "MPI_SIZET_DATA_TYPE is undefined"
#endif

// bitset to store code for space filling curves
static constexpr DefUint kSFBitsetBit = 64;
using DefSFBitset = std::bitset<kSFBitsetBit>;
class HashFunc {
 public:
    DefSFCodeToUint operator()(const DefSFBitset&bit_set_in) const noexcept {
        DefSFCodeToUint hashVal = bit_set_in.to_ullong();
        return hashVal;
    }
};
template <typename value>
using DefMap = std::unordered_map<DefSFBitset, value, HashFunc>;

// numerical parameters
const DefReal kEps = 1.e-10;
/**< Constant: a small number for comparing floats*/
const DefReal kPi = 4. * atan(1.);  /**< Constant: pi*/

// globle constants
static constexpr DefUint kXIndex = 0, kYIndex = 1, kZIndex = 2;
/**< indices of vectors for x, y, z diretions*/
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_DEF_LIB_H_
