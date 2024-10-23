#!/bin/sh

# Test test inter mtt coding.

set -eu

. "${0%/*}/util.sh"

common_args='264x130 10 yuv420p --preset=ultrafast --threads=0 --no-cpuid --no-wpp --fast-residual-cost 0'

valgrind_test $common_args --rd=0 --mtt-depth-inter 1 --pu-depth-inter 2-3
valgrind_test $common_args --rd=3 --mtt-depth-inter 1 --pu-depth-inter 0-5
valgrind_test $common_args --rd=3 --mtt-depth-inter 3 --pu-depth-inter 0-8
valgrind_test $common_args --rd=3 --mtt-depth-inter 3 --pu-depth-inter 0-8 --mts inter
