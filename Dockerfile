# A simple Dockerfile for building uvg266 from the git repository
# Example build command when in this directory: docker build -t uvg266 .
#
# Example usage
# Run with an input YUV file and output HEVC binary file
#     docker run -i -a STDIN -a STDOUT uvg266 -i - --input-res=320x240 -o - < testfile_320x240.yuv > out.265
#
# Use libav or ffmpeg to input (almost) any format and convert it to YUV420 for uvg266, audio is disabled
#
#     RESOLUTION=`avconv -i input.avi 2>&1 | grep Stream | grep -oP ', \K[0-9]+x[0-9]+'`
#     avconv -i input.avi -an -f rawvideo -pix_fmt yuv420p - | docker run -i -a STDIN -a STDOUT uvg266 -i - --wpp --threads=8 --input-res=$RESOLUTION --preset=ultrafast -o - > output.266
#  or
#     RESOLUTION=`ffmpeg -i input.avi 2>&1 | grep Stream | grep -oP ', \K[0-9]+x[0-9]+'`
#     ffmpeg -i input.avi -an -f rawvideo -pix_fmt yuv420p - | docker run -i -a STDIN -a STDOUT uvg266 -i - --wpp --threads=8 --input-res=$RESOLUTION --preset=ultrafast -o - > output.266
#

# Use Ubuntu 22.04 as a base for now
FROM ubuntu:22.04

MAINTAINER Marko Viitanen <fador@iki.fi>

# List of needed packages to be able to build uvg266 with autotools
ENV REQUIRED_PACKAGES automake autoconf libtool m4 build-essential git cmake pkgconf

COPY . uvg266
# Run all the commands in one RUN so we don't have any extra history
# data in the image.
RUN apt-get update \
    && apt-get install -y $REQUIRED_PACKAGES \
    && cd uvg266/build \
    && cmake -DUSE_SHARED_LIB=OFF .. \
    && make\
    && make install \
    && AUTOINSTALLED_PACKAGES=`apt-mark showauto` \
    && apt-get remove --purge --force-yes -y $REQUIRED_PACKAGES $AUTOINSTALLED_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

ENTRYPOINT ["uvg266"]
CMD ["--help"]
