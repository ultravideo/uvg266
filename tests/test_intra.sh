#!/bin/sh

# Test all-intra coding.

set -eu

. "${0%/*}/util.sh"

common_args='256x128 10 -p1 --preset=ultrafast --threads=0 --no-wpp --no-tmvp --no-deblock --sao=0 --cpuid=0'
valgrind_test $common_args --rd=1
valgrind_test $common_args --rd=2 --no-transform-skip --qp 37
