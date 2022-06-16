#!/bin/sh

LANG=C
set -e

cd "$(dirname "$0")"

date="$(date +"%B %Y")"
version="$(awk 'match($0,/VERSION [0-9]\.[0-9]\.[0-9]/) {print $2}' ../CMakeLists.txt)"
manpage_file=../doc/uvg266.1

cat <<EOF> $manpage_file
.TH UVG266 "1" "$date" "uvg266 v$version" "User Commands"
.SH NAME
uvg266 \- open source VVC encoder
.SH SYNOPSIS
\fBuvg266 \fR\-i <input> \-\-input\-res <width>x<height> \-o <output>
.SH DESCRIPTION
EOF

../bin/uvg266 --help 2>&1 | tail -n+5 | \
  sed 's| : |\n|g;
       s| :$||g;
       s|^      --|.TP\n\\fB--|g;
       s|^  -|.TP\n\\fB-|g;
       s|^                               ||g;
       s|-|\\-|g;
       s|, \\-\\-|\\fR, \\fB\\-\\-|g;' \
  >> $manpage_file

for s in Required Presets Input Options "Video structure" "Compression tools" "Parallel processing" "Video Usability Information"; do
  sed -i "s|^${s}:|.SS \"${s}:\"|g" $manpage_file
done
