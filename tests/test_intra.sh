#!/bin/sh

# Test all-intra coding.

set -eu

. "${0%/*}/util.sh"

common_args='264x130 10 -p1 --threads=0 --no-wpp --no-rdoq --no-deblock --no-sao --no-signhide --cpuid=0'
valgrind_test $common_args --rd=1
valgrind_test $common_args --rd=2 --no-transform-skip
