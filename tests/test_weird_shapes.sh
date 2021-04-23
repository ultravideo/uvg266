#!/bin/sh

set -eu
. "${0%/*}/util.sh"

valgrind_test  16x16  10 yuv420p --threads=0 --no-wpp --preset=veryslow
valgrind_test 256x16  10 yuv420p --threads=0 --no-wpp --preset=veryslow
valgrind_test  16x256 10 yuv420p --threads=0 --no-wpp --preset=veryslow
