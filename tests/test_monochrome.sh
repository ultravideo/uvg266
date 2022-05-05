#!/bin/sh

# Test all-intra coding.

set -eu

. "${0%/*}/util.sh"

common_args='256x128 10 gray --hash none --preset=ultrafast --input-format P400 --threads=0 --no-wpp --no-tmvp --no-deblock --sao=0 --cpuid=0 --pu-depth-intra 0-4'
valgrind_test $common_args --rd=1 -p1 
valgrind_test $common_args --rd=2 -p1 
valgrind_test $common_args --rd=2 -p1 --deblock 0:0
