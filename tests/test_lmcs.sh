#!/bin/sh

# Test lmcs

set -eu

. "${0%/*}/util.sh"

common_args='256x128 10 yuv420p -p1 --preset=ultrafast --threads=0 --no-wpp --no-tmvp --no-deblock --sao=0 --cpuid=0 --pu-depth-intra 0-3 --lmcs'
valgrind_test $common_args --deblock --sao=full --alf=full

