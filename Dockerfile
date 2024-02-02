# Set the base image to be the latest Ubuntu image
FROM ubuntu:20.04

# Set the working directory to be build
WORKDIR /build

# Install dependencies
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata &&\
    apt-get install -y zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    git \
    cmake \
    python3 \
    curl \
    gcc-9 \
    g++-9 \
    bzip2 \
    autoconf \
    automake \ 
    && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9

# Install Moni

ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/build/moni-align/build/thirdparty/lib"

RUN git clone https://github.com/maxrossi91/moni-align.git &&\
    cd moni-align &&\
    git checkout develop &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    make || true &&\
    cmake .. &&\
    make 

WORKDIR /mnt
ENV PATH /build/moni-align/build:$PATH