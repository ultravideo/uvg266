#pragma once

#include <functional>
#include <map>
#include <string>

struct Benchmark
{
  std::string name;
  std::function<int(int, char**)> init;
};

inline auto get_benchmarks() -> std::map<std::string, Benchmark>&
{
  static std::map<std::string, Benchmark> instance;
  return instance;
}

inline auto register_benchmark(Benchmark&& benchmark) -> bool
{
  get_benchmarks()[benchmark.name] = benchmark;
  return true;
}

#define UVG_BENCHMARK_CONCAT_INNER(a, b) a##b
#define UVG_BENCHMARK_CONCAT(a, b) UVG_BENCHMARK_CONCAT_INNER(a, b)

#if defined(__COUNTER__)
  #define UVG_BENCHMARK_UNIQUE_NAME(base) UVG_BENCHMARK_CONCAT(base, __COUNTER__)
#else
  #define UVG_BENCHMARK_UNIQUE_NAME(base) UVG_BENCHMARK_CONCAT(base, __LINE__)
#endif

#define UVG_BENCHMARK_REGISTER(name, func) \
  static const bool UVG_BENCHMARK_UNIQUE_NAME(_auto_reg_) = [] { \
    return register_benchmark({(name), (func)}); \
  }();
