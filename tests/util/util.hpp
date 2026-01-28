#pragma once

extern "C" {
#include "global.h"
}

#include <array>
#include <cassert>
#include <ctime>
#include <limits>
#include <random>
#include <vector>

#ifndef UVG_NO_TESTS
#include <catch2/catch_get_random_seed.hpp>
#endif

namespace util {

inline auto rng() -> std::mt19937&
{
#ifndef UVG_NO_TESTS
  static const auto seed = static_cast<std::mt19937::result_type>(Catch::getSeed());
#else
  static const auto seed = static_cast<std::mt19937::result_type>(time(nullptr));
#endif
  static std::mt19937 instance(seed);
  return instance;
}

template<typename T>
auto random_int() -> T
{
  const auto min = std::numeric_limits<T>::min();
  const auto max = std::numeric_limits<T>::max();
  return std::uniform_int_distribution<T>(min, max)(rng());
}

template<typename T>
auto random_nonnegative() -> T
{
  const auto max = std::numeric_limits<T>::max();
  return std::uniform_int_distribution<T>(0, max)(rng());
}

template<typename T>
auto random_positive() -> T
{
  const auto max = std::numeric_limits<T>::max();
  return std::uniform_int_distribution<T>(1, max)(rng());
}

template<typename T>
auto random_bounded(T min, T max) -> T
{
  T value = min + (random_positive<T>() % (max - min + 1));
  return value;
}

// Specialization to satisfy MSVC.
template<>
inline auto random_int<uint8_t>() -> uint8_t
{
  return static_cast<uint8_t>(random_int<uint64_t>());
}

// Specialization to satisfy MSVC.
template<>
inline auto random_int<int8_t>() -> int8_t
{
  return static_cast<int8_t>(random_int<int64_t>());
}

template<typename T, size_t S>
auto random_int_array() -> std::array<T, S>
{
  static_assert(std::is_integral_v<T>);

  std::array<T, S> array;
  for (size_t i = 0; i < S; ++i) {
    array[i] = random_int<T>();
  }
  return array;
}

template<typename T>
auto random_int_vector(size_t size) -> std::vector<T>
{
  static_assert(std::is_integral_v<T>);

  std::vector<T> vector;
  vector.resize(size);
  for (size_t i = 0; i < size; ++i) {
    vector[i] = random_int<T>();
  }
  return vector;
}

template<typename T>
auto random_of(std::initializer_list<T> choices) -> T
{
    assert(choices.size() > 0);
    return *(choices.begin() + random_bounded<size_t>(0, choices.size() - 1));
}

template<typename T>
auto random_of(const std::vector<T>& choices) -> T
{
    return choices[random_bounded<size_t>(0, choices.size() - 1)];
}

template<typename T, size_t N>
void fill_sequential(std::array<T, N>& array)
{
  for (size_t i = 0; i < N; ++i) {
    array[i] = static_cast<T>(i);
  }
}

template<typename T, size_t N>
void fill_random(std::array<T, N>& array)
{
  for (size_t i = 0; i < N; ++i) {
    array[i] = random_int<T>();
  }
}

template<typename T, size_t N>
void fill_random_bounded(std::array<T, N>& array, T min, T max)
{
  for (size_t i = 0; i < N; ++i) {
    array[i] = random_bounded<T>(min, max);
  }
}

template<typename T, size_t N>
void fill_constant(std::array<T, N>& array, typename std::array<T, N>::value_type value)
{
  for (size_t i = 0; i < N; ++i) {
    array[i] = value;
  }
}

}
