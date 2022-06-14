#!/bin/sh

# Test RDOQ, SAO, deblock and signhide and subme.

set -eu
. "${0%/*}/util.sh"

common_args='264x130 10 yuv420p -p0 -r1 --threads=2 --wpp --owf=1 --rd=0 --pu-depth-inter 0-3 --no-bipred --no-tmvp --gop=0'

valgrind_test $common_args --no-rdoq --no-deblock --no-sao --no-signhide --subme=1 --pu-depth-intra=2-3
valgrind_test $common_args --no-rdoq --no-signhide --subme=0 --bipred
valgrind_test $common_args --rdoq --no-deblock --no-sao --subme=0
valgrind_test $common_args --gop=8 --subme=4 --bipred --tmvp
valgrind_test $common_args --transform-skip --tr-skip-max-size=5
valgrind_test $common_args --vaq=8
valgrind_test $common_args --vaq=8 --bitrate 350000
valgrind_test $common_args --vaq=8 --rc-algorithm oba --bitrate 350000