#!/bin/sh

# Test RDOQ, SAO, deblock and signhide and subme.

set -eu
. "${0%/*}/util.sh"

common_args='264x130 10 yuv420p -p0 -r1 --threads=2 --wpp --owf=1 --rd=0'

valgrind_test $common_args --no-rdoq --no-deblock --no-sao --no-signhide --subme=1 --pu-depth-intra=2-3
valgrind_test $common_args --no-rdoq --no-signhide --subme=0
valgrind_test $common_args --rdoq --no-deblock --no-sao --subme=0
valgrind_test $common_args --vaq=8
valgrind_test $common_args --vaq=8 --bitrate 3500
valgrind_test $common_args --vaq=8 --rc-algorithm oba --bitrate 3500
valgrind_test $common_args --lmcs
valgrind_test $common_args --jccr
valgrind_test $common_args --jccr --rdoq --rd=2 --mts=full
