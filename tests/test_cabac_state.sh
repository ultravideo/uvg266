#!/bin/sh

set -eu

. "${0%/*}/util.sh"

cabacfile="$(mktemp)"

valgrind_test 256x128 10 yuv420p --preset veryslow --rd 3 --mip --jccr --mrl --lfnst -p 1 --owf 0 --no-wpp --cabac-debug-file="${cabacfile}"
python3 check_cabac_state_consistency.py "${cabacfile}"

valgrind_test 256x128 10 yuv420p --preset veryslow --rd 3 --mip --jccr --mrl --lfnst --dual-tree -p 1 --owf 0 --no-wpp --cabac-debug-file="${cabacfile}"
python3 check_cabac_state_consistency.py "${cabacfile}"

rm -rf "${cabacfile}"

