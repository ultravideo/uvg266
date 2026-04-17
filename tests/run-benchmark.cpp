#include "run-benchmark.hpp"

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <benchmark/benchmark.h>
#include <CLI/CLI.hpp>

static auto print_usage(std::string_view program_name, std::ostream& to) -> void
{
  to << "Usage:\n  " << program_name << " <benchmark> [options] [-- Google Benchmark options]\n";
  to << '\n';
  to << "  where benchmark is one of:\n";
  for (const auto& [name, benchmark] : get_benchmarks()) {
    to << "    " << name << '\n';
  }
}

auto main(int argc, char* argv[]) -> int
{
  benchmark::MaybeReenterWithoutASLR(argc, argv);

  std::vector<char*> our_args = {argv[0]};
  std::vector<char*> their_args = {argv[0]};

  const std::string_view arg_separator = "--";
  bool arg_separator_found = false;
  for (int i = 1; i < argc; ++i) {
    if (argv[i] == arg_separator) {
      arg_separator_found = true;
      continue;
    }
    if (!arg_separator_found) {
      our_args.push_back(argv[i]);
    } else {
      their_args.push_back(argv[i]);
    }
  }

  int their_argc = static_cast<int>(their_args.size());
  char** their_argv = their_args.data();
  benchmark::Initialize(&their_argc, their_argv);
  benchmark::ReportUnrecognizedArguments(their_argc, their_argv);

  int our_argc = static_cast<int>(our_args.size());
  char** our_argv = our_args.data();

  if (our_argc < 2) {
    print_usage(our_argv[0], std::cerr);
    return EXIT_FAILURE;
  }

  const std::string first_arg = our_argv[1];
  if (first_arg == "--help" || first_arg == "-h") {
    print_usage(our_argv[0], std::cout);
    return EXIT_SUCCESS;
  }

  const auto& benchmarks = get_benchmarks();
  const auto iter = benchmarks.find(first_arg);
  if (iter == benchmarks.end()) {
    std::cerr << "No such benchmark: " << first_arg << '\n';
    print_usage(our_argv[0], std::cerr);
    return EXIT_FAILURE;
  }

  const auto& [name, benchmark] = *iter;
  benchmark.init(our_argc - 1, our_argv + 1);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return EXIT_SUCCESS;
}
