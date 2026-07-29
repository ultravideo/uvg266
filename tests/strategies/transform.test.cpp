#include "cu.h"
extern "C" {
#include "strategies/strategies-dct.h"
#include "strategies/generic/dct-generic.h"
#include "uvg266.h"
}

#include "../util/util.hpp"

#include <array>
#include <thread>

// ---------------------------------------------------------
// TESTS
// ---------------------------------------------------------

#ifndef UVG_NO_TESTS

#include <catch2/catch_test_macros.hpp>

namespace {

struct Parameters
{
  int8_t width = 0;
  int8_t height = 0;
  int8_t bitdepth = 0;
  int8_t color = 0;
  int8_t mts_type = 0;
  uint8_t tr_idx = 0;
  int8_t isp_mode = 0;
};

}

static auto test_parameters() -> std::vector<Parameters>
{
  std::vector<Parameters> parameters;

  for (int8_t width : {1, 2, 4, 8, 16, 32}) {
    for (int8_t height : {1, 2, 4, 8, 16, 32}) {
      if ((width == 1 && height < 16) || (height == 1 && width < 16) || (width == 2 && height == 2)) {
        continue;
      }
      for (int8_t bitdepth : {8, 10}) {
        for (int8_t color : {COLOR_Y, COLOR_U}) {
          for (int8_t mts_type : {UVG_MTS_IMPLICIT, UVG_MTS_BOTH}) {
            std::vector<uint8_t> tr_idxs = {MTS_DCT2_DCT2};
            if ((width >= 4 && height >= 4) && mts_type == UVG_MTS_BOTH && color == COLOR_Y) {
              tr_idxs = {
                MTS_DCT2_DCT2,
                MTS_DST7_DST7,
                MTS_DCT8_DST7,
                MTS_DST7_DCT8,
                MTS_DCT8_DCT8,
              };
            }
            const int8_t isp_mode = (color == COLOR_Y) && ((width < 4) || (height < 4));
            for (uint8_t tr_index : tr_idxs) {
              parameters.push_back({
                width,
                height,
                bitdepth,
                color,
                mts_type,
                tr_index,
                isp_mode
              });
            }
          }
        }
      }
    }
  }

  return parameters;
}

static void run_test_mts_dct(mts_dct_func* tested_impl, mts_dct_func* golden_impl)
{
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> input;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> actual;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> expected;

  for (const auto& [width, height, bitdepth, color, mts_type, tr_index, isp_mode] : test_parameters()) {
    INFO("width=" << static_cast<int>(width)
      << " height=" << static_cast<int>(height)
      << " bitdepth=" << static_cast<int>(bitdepth)
      << " color=" << static_cast<int>(color)
      << " mts_type=" << static_cast<int>(mts_type)
      << " tr_index=" << static_cast<int>(tr_index));
    for (int i = 0; i < 5; ++i) {
      const int16_t min_input = bitdepth == 8 ? -255 : -1023;
      const int16_t max_input = bitdepth == 8 ? 255 : 1023;
      util::fill_random_bounded<int16_t>(input, min_input, max_input);
      util::fill_constant(actual, -1);
      util::fill_constant(expected, -1);
      cu_info_t tu = {};
      tu.type = CU_INTRA;
      tu.tr_idx = tr_index;
      tu.intra.isp_mode = isp_mode;
      tested_impl(bitdepth, (color_t)color, &tu, width, height, input.data(), actual.data(), mts_type);
      golden_impl(bitdepth, (color_t)color, &tu, width, height, input.data(), expected.data(), mts_type);
      REQUIRE(actual == expected);
    }
  }
}

static void run_test_mts_idct(mts_dct_func* tested_impl, mts_dct_func* golden_impl)
{
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> input;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> actual;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> expected;

  for (const auto& [width, height, bitdepth, color, mts_type, tr_index, isp_mode] : test_parameters()) {
    INFO("width=" << static_cast<int>(width)
      << " height=" << static_cast<int>(height)
      << " bitdepth=" << static_cast<int>(bitdepth)
      << " color=" << static_cast<int>(color)
      << " mts_type=" << static_cast<int>(mts_type)
      << " tr_index=" << static_cast<int>(tr_index));
    for (int i = 0; i < 5; ++i) {
      const int16_t min_input = INT16_MIN;
      const int16_t max_input = INT16_MAX;
      util::fill_random_bounded<int16_t>(input, min_input, max_input);
      util::fill_constant(actual, -1);
      util::fill_constant(expected, -1);
      cu_info_t tu = {};
      tu.type = CU_INTRA;
      tu.tr_idx = tr_index;
      tu.intra.isp_mode = isp_mode;

      if (tr_index != MTS_DCT2_DCT2) {
        for (int y = 0; y < height; ++y) {
          for (int x = 0; x < width; ++x) {
            if (y >= 16 || x >= 16) {
              input[y * width + x] = 0;
            }
          }
        }
      }

      tested_impl(bitdepth, (color_t)color, &tu, width, height, input.data(), actual.data(), mts_type);
      golden_impl(bitdepth, (color_t)color, &tu, width, height, input.data(), expected.data(), mts_type);
      REQUIRE(actual == expected);
    }
  }
}

static void run_test_dct(dct_func* tested_impl, dct_func* golden_impl)
{
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> input;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> actual;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> expected;

  for (int bitdepth : {8, 10}) {
    INFO("bitdepth=" << bitdepth);
    const int16_t min_input = bitdepth == 8 ? -255 : -1023;
    const int16_t max_input = bitdepth == 8 ? 255 : 1023;
    util::fill_random_bounded<int16_t>(input, min_input, max_input);
    util::fill_constant(actual, -1);
    util::fill_constant(expected, -1);
    tested_impl(static_cast<int8_t>(bitdepth), input.data(), actual.data());
    golden_impl(static_cast<int8_t>(bitdepth), input.data(), expected.data());
    REQUIRE(actual == expected);
  }
}

static void run_test_idct(dct_func* tested_impl, dct_func* golden_impl)
{
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> input;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> actual;
  ALIGNED(SIMD_ALIGNMENT) std::array<int16_t, TR_MAX_WIDTH * TR_MAX_WIDTH> expected;

  for (int bitdepth : {8, 10}) {
    INFO("bitdepth=" << bitdepth);
    const int16_t min_input = INT16_MIN;
    const int16_t max_input = INT16_MAX;
    util::fill_random_bounded<int16_t>(input, min_input, max_input);
    util::fill_constant(actual, -1);
    util::fill_constant(expected, -1);
    tested_impl(static_cast<int8_t>(bitdepth), input.data(), actual.data());
    golden_impl(static_cast<int8_t>(bitdepth), input.data(), expected.data());
    REQUIRE(actual == expected);
  }
}

#endif // UVG_NO_TESTS

// ---------------------------------------------------------
// BENCHMARKS
// ---------------------------------------------------------

#ifndef UVG_NO_BENCHMARKS

#include "../run-benchmark.hpp"

#include <benchmark/benchmark.h>
#include <CLI/CLI.hpp>

template<mts_idx tr_idx, auto mts_dct_func>
static void benchmark_run(benchmark::State& state)
{
  ALIGNED(SIMD_ALIGNMENT) int16_t input[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(SIMD_ALIGNMENT) int16_t output[TR_MAX_WIDTH * TR_MAX_WIDTH];

  const auto bitdepth = static_cast<int8_t>(state.range(0));
  const auto width = static_cast<int8_t>(state.range(1));
  const auto height = static_cast<int8_t>(state.range(2));

  cu_info_t tu = {};
  tu.type = CU_INTRA;
  tu.tr_idx = tr_idx;

  for (auto _ : state) {
    mts_dct_func(bitdepth, COLOR_Y, &tu, width, height, input, output, UVG_MTS_BOTH);
  }
}

template<mts_idx tr_idx, auto mts_dct_func>
static auto benchmark_init(int argc, char** argv) -> int
{
  std::string name = argv[0];

  struct BenchmarkArgs {
    int bitdepth = 8;
    std::vector<std::pair<size_t, size_t>> size = {
      { 1, 16}, { 1, 32},
      { 2,  4}, { 2,  8}, { 2, 16}, { 2, 32},
      { 4,  2}, { 4,  4}, { 4,  8}, { 4, 16}, { 4, 32},
      { 8,  2}, { 8,  4}, { 8,  8}, { 8, 16}, { 8, 32},
      {16,  1}, {16,  2}, {16,  4}, {16,  8}, {16, 16}, {16, 32},
      {32,  1}, {32,  2}, {32,  4}, {32,  8}, {32, 16}, {32, 32},
    };
    unsigned threads = std::thread::hardware_concurrency();
  } args;

  CLI::App cli;
  cli.add_option("-b,--bitdepth", args.bitdepth);
  cli.add_option("-s,--size", args.size);
  cli.add_option("-t,--threads", args.threads);

  CLI11_PARSE(cli, argc, argv);

  auto* benchmark = benchmark::RegisterBenchmark(name, benchmark_run<tr_idx, mts_dct_func>);
  benchmark->ArgNames({"bitdepth", "width", "height"});

  for (const auto& [width, height] : args.size) {
    if ((width >= 4 && height >= 4) || (tr_idx == MTS_DCT2_DCT2)) {
      benchmark->Args({static_cast<int64_t>(args.bitdepth), static_cast<int64_t>(width), static_cast<int64_t>(height)});
    }
  }

  if (args.threads > 1) {
    benchmark->Threads(static_cast<int>(args.threads));
  }

  return EXIT_SUCCESS;
}

UVG_BENCHMARK_REGISTER("dct2_forward_generic", (benchmark_init<MTS_DCT2_DCT2, uvg_mts_dct_generic>));
UVG_BENCHMARK_REGISTER("dct2_inverse_generic", (benchmark_init<MTS_DCT2_DCT2, uvg_mts_idct_generic>));
UVG_BENCHMARK_REGISTER("dct8_forward_generic", (benchmark_init<MTS_DCT8_DCT8, uvg_mts_dct_generic>));
UVG_BENCHMARK_REGISTER("dct8_inverse_generic", (benchmark_init<MTS_DCT8_DCT8, uvg_mts_idct_generic>));
UVG_BENCHMARK_REGISTER("dst7_forward_generic", (benchmark_init<MTS_DST7_DST7, uvg_mts_dct_generic>));
UVG_BENCHMARK_REGISTER("dst7_inverse_generic", (benchmark_init<MTS_DST7_DST7, uvg_mts_idct_generic>));

#endif // UVG_NO_BENCHMARKS
