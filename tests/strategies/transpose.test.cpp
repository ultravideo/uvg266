extern "C" {
#include "build_config.h"
#include "strategies/generic/transpose-generic.h"
#include "strategies/neon/transpose-neon.h"
#include "strategies/strategies-transpose.h"
}

#include "../util/util.hpp"

#include <array>
#include <thread>

namespace {

constexpr int MAX_WIDTH = 64;
constexpr int MAX_HEIGHT = 64;

struct Parameters
{
  int width = 0;
  int height = 0;
};

template<typename T>
inline void transpose_reference(const T* src, T* dst, int width, int height)
{
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      dst[y + x * height] = src[x + y * width];
    }
  }
}

void transpose_8bit_reference(const void* src, void* dst, int width, int height)
{
  return transpose_reference(static_cast<const uint8_t*>(src), static_cast<uint8_t*>(dst), width, height);
}

void transpose_16bit_reference(const void* src, void* dst, int width, int height)
{
  return transpose_reference(static_cast<const uint16_t*>(src), static_cast<uint16_t*>(dst), width, height);
}

void transpose_px_reference(const uvg_pixel* src, uvg_pixel* dst, int width, int height)
{
  return transpose_reference(src, dst, width, height);
}

} // namespace

// ---------------------------------------------------------
// TESTS
// ---------------------------------------------------------

#ifndef UVG_NO_TESTS

#include <catch2/catch_test_macros.hpp>

namespace {

auto test_parameters() -> std::vector<Parameters>
{
  std::vector<Parameters> parameters;

  // Test powers of two.
  for (int width = 1; width <= MAX_WIDTH; width *= 2) {
    for (int height = 1; height <= MAX_HEIGHT; height *= 2) {
      parameters.push_back({width, height});
    }
  }

  // Test a bunch of random sizes.
  for (int i = 0; i < 100; ++i) {
    int width = util::random_bounded<int>(1, MAX_WIDTH);
    int height = util::random_bounded<int>(1, MAX_HEIGHT);
    parameters.push_back({width, height});
  }

  return parameters;
}

template<typename T, typename func>
void run_test(func* tested_impl, func* golden_impl)
{
  ALIGNED(SIMD_ALIGNMENT) std::array<T, MAX_WIDTH * MAX_HEIGHT> input;
  ALIGNED(SIMD_ALIGNMENT) std::array<T, MAX_WIDTH * MAX_HEIGHT> actual;
  ALIGNED(SIMD_ALIGNMENT) std::array<T, MAX_WIDTH * MAX_HEIGHT> expected;
  util::fill_sequential(input);
  util::fill_constant(actual, 0);
  util::fill_constant(expected, 0);

  for (const auto& [width, height] : test_parameters()) {
    INFO("size=" << width << "x" << height);
    tested_impl(actual.data(), input.data(), width, height);
    golden_impl(expected.data(), input.data(), width, height);
    REQUIRE(memcmp(actual.data(), expected.data(), static_cast<size_t>(width * height) * sizeof(T)) == 0);
  }
}

} // namespace

TEST_CASE("uvg_transpose_8bit_generic")
{
  run_test<uint8_t, transpose_8bit_func>(uvg_transpose_8bit_generic, transpose_8bit_reference);
}

TEST_CASE("uvg_transpose_16bit_generic")
{
  run_test<uint16_t, transpose_16bit_func>(uvg_transpose_16bit_generic, transpose_16bit_reference);
}

TEST_CASE("uvg_transpose_px_generic")
{
  run_test<uvg_pixel, transpose_px_func>(uvg_transpose_px_generic, transpose_px_reference);
}

#if UVG_HAVE_ARM_NEON

TEST_CASE("uvg_transpose_8bit_neon")
{
  run_test<uint8_t, transpose_8bit_func>(uvg_transpose_8bit_neon, transpose_8bit_reference);
}

TEST_CASE("uvg_transpose_16bit_neon")
{
  run_test<uint16_t, transpose_16bit_func>(uvg_transpose_16bit_neon, transpose_16bit_reference);
}

TEST_CASE("uvg_transpose_px_neon")
{
  run_test<uvg_pixel, transpose_px_func>(uvg_transpose_px_neon, transpose_px_reference);
}

#endif

#endif // UVG_NO_TESTS

// ---------------------------------------------------------
// BENCHMARKS
// ---------------------------------------------------------

#ifndef UVG_NO_BENCHMARKS

#include "../run-benchmark.hpp"

#include <benchmark/benchmark.h>
#include <CLI/CLI.hpp>

namespace {

template<typename elem_type, auto transpose_func>
auto benchmark_run(benchmark::State& state) -> void
{
  ALIGNED(SIMD_ALIGNMENT) std::array<elem_type, MAX_WIDTH * MAX_HEIGHT> in;
  ALIGNED(SIMD_ALIGNMENT) std::array<elem_type, MAX_WIDTH * MAX_HEIGHT> out;
  util::fill_sequential(in);
  util::fill_constant(out, 0);

  const int width = static_cast<int>(state.range(0));
  const int height = static_cast<int>(state.range(1));

  for (auto _ : state) {
    transpose_func(out.data(), in.data(), width, height);
  }
}

template<typename elem_type, auto transpose_func>
auto benchmark_init(int argc, char** argv) -> int
{
  std::string name = argv[0];

  static constexpr auto default_sizes = [](){
    std::vector<std::pair<unsigned, unsigned>> sizes;
    for (unsigned w = 1; w <= MAX_WIDTH; w *= 2) {
      for (unsigned h = 1; h <= MAX_WIDTH; h *= 2) {
        sizes.push_back({w, h});
      }
    }
    return sizes;
  };

  struct BenchmarkArgs {
    std::vector<std::pair<unsigned, unsigned>> size = default_sizes();
    unsigned threads = std::thread::hardware_concurrency();
  } args;

  CLI::App cli;
  cli.add_option("-s,--size", args.size);
  cli.add_option("-t,--threads", args.threads);

  CLI11_PARSE(cli, argc, argv);

  auto* benchmark = benchmark::RegisterBenchmark(name, benchmark_run<elem_type, transpose_func>);
  benchmark->ArgNames({"width", "height"});

  for (const auto& [width, height] : args.size) {
    benchmark->Args({static_cast<int64_t>(width), static_cast<int64_t>(height)});
  }

  if (args.threads > 1) {
    benchmark->Threads(static_cast<int>(args.threads));
  }

  return EXIT_SUCCESS;
}

} // namespace

UVG_BENCHMARK_REGISTER("uvg_transpose_8bit_generic", (benchmark_init<uint8_t, uvg_transpose_8bit_generic>));
UVG_BENCHMARK_REGISTER("uvg_transpose_16bit_generic", (benchmark_init<uint16_t, uvg_transpose_16bit_generic>));

#if UVG_HAVE_ARM_NEON
UVG_BENCHMARK_REGISTER("uvg_transpose_8bit_neon", (benchmark_init<uint8_t, uvg_transpose_8bit_neon>));
UVG_BENCHMARK_REGISTER("uvg_transpose_16bit_neon", (benchmark_init<uint16_t, uvg_transpose_16bit_neon>));
#endif

#endif // UVG_NO_BENCHMARKS
