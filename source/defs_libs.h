//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file defs_libs.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief contain definations and libraries will be used in all modules.
*/

#ifndef SOURCE_DEFS_LIBS_H_
#define SOURCE_DEFS_LIBS_H_
// libraries
#include <bitset>
#include <cstdint>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "./config.h"  // configuration file generated
namespace rootproject {
using DefReal = double_t;
using DefInt = int32_t;
using DefUint = uint32_t;
using DefSizet = size_t;
// unsigned integer for indices will be long like indices of nodes and vertices
using DefAmrLUint = uint32_t;
// unsigned integer for space filling code
using DefSFCodeToUint = uint64_t;


#ifdef ENABLE_MPI
// data type for mpi communication
#define MPI_REAL_DATA_TYPE MPI_DOUBLE
#define MPI_INT_DATA_TYPE MPI_INT32_T
#define MPI_UINT_DATA_TYPE MPI_INT32_T
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
#define MPI_AMR_LUINT_TYPE MPI_INT32_T
#define MPI_CODE_UINT_TYPE MPI_UINT64_T

#endif

// bitset to store code for space filling curves
static constexpr DefInt kSFBitsetBit = 64;
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
const DefReal kEps = DefReal(1.e-6);
/**< Constant: a small number for comparing floats*/
const DefReal kPi = DefReal(4.) * atan(DefReal(1.));  /**< Constant: pi*/

// global constants
static constexpr DefInt kXIndex = 0, kYIndex = 1, kZIndex = 2;
/**< indices of vectors for x, y, z directions*/
}  // end namespace rootproject
#endif  // SOURCE_DEFS_LIBS_H_
