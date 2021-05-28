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

wget http://ultravideo.fi/ubuntu-vtm-13.0.tgz
sha256sum -c - << EOF
1751bba04611791040a978015da8bb537fcd89c89453c0cdc1ec2719aac2a1a8  ubuntu-vtm-13.0.tgz
EOF
tar xf ubuntu-vtm-13.0.tgz
cp DecoderApp "${HOME}/bin/DecoderAppStatic"
chmod +x "${HOME}/bin/DecoderAppStatic"

export PATH="${HOME}/bin:${PATH}"
