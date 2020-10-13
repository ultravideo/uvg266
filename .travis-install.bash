#!/bin/bash

# Download FFmpeg and HM decoder and place them in $PATH.

set -euvo pipefail

mkdir -p "${HOME}/bin"

wget http://ultravideo.fi/ffmpeg-release-4.2.1-32bit-static.tar.xz
sha256sum -c - << EOF
226f55f8a94d71f3d231a20fe59fcbb7f6100cabf663f9bcb887d17b332a91c5  ffmpeg-release-4.2.1-32bit-static.tar.xz
EOF
tar xf ffmpeg-release-4.2.1-32bit-static.tar.xz
cp ffmpeg-4.2.1-i686-static/ffmpeg "${HOME}/bin/ffmpeg"
chmod +x "${HOME}/bin/ffmpeg"

wget http://ultravideo.fi/ubuntu-vtm-ea11cf4b.tgz
sha256sum -c - << EOF
3216516a397d4a7ed3eb18279c53e69ee766fe5cdda7c6f3b2f0b0de951c716c  ubuntu-vtm-ea11cf4b.tgz
EOF
tar xf ubuntu-vtm-ea11cf4b.tgz
cp DecoderAppStatic "${HOME}/bin/DecoderAppStatic"
chmod +x "${HOME}/bin/DecoderAppStatic"

export PATH="${HOME}/bin:${PATH}"
