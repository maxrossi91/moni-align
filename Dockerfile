# Set the base image to be the latest Ubuntu image
FROM ubuntu:latest

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

# Download and install the main branch of the sdsl-lite library
RUN git clone https://github.com/simongog/sdsl-lite.git &&\
    cd sdsl-lite &&\
    ./install.sh

# Install htslib/1.15 library
RUN curl -LJO https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 &&\
    tar -xvjf htslib-1.15.tar.bz2 &&\
    cd htslib-1.15 &&\
    autoreconf -i &&\
    ./configure &&\
    make &&\ 
    make install &&\
    rm ../htslib-1.15.tar.bz2

# Running ldconfig after install so newly installed libraries are found
RUN ldconfig

# Install Moni
RUN git clone https://github.com/maxrossi91/moni-align.git &&\
    cd moni-align &&\
    git checkout develop &&\
    mkdir build &&\
    cd build &&\
    cmake -DCMAKE_LIBRARY_PATH="/build/htslib-1.15/;/root/lib/" -DCMAKE_INCLUDE_PATH="/build/htslib-1.15/;/root/include/" -DCMAKE_PREFIX_PATH="`realpath thirdparty`" .. || true &&\
    cd _deps/pfp-build &&\
    cmake -DCMAKE_LIBRARY_PATH="/build/htslib-1.15/;/root/lib/" -DCMAKE_INCLUDE_PATH="/build/htslib-1.15/;/root/include/" -DCMAKE_PREFIX_PATH="`realpath thirdparty`" -DCMAKE_INSTALL_PREFIX="`realpath ../../thirdparty`" ../pfp-src &&\
    make &&\
    make install &&\
    cd ../../ &&\
    cmake -DCMAKE_LIBRARY_PATH="/build/htslib-1.15/;/root/lib/" -DCMAKE_INCLUDE_PATH="/build/htslib-1.15/;/root/include/" -DCMAKE_PREFIX_PATH="`realpath thirdparty`" .. &&\
    make

# To ensure paths in Moni python script are set correct
ENV PATH="$PATH:/build/moni-align/build/thirdparty/bin"
RUN cd moni-align &&\
    cd build &&\
    cmake -DCMAKE_LIBRARY_PATH="/build/htslib-1.15/;/root/lib/" -DCMAKE_INCLUDE_PATH="/build/htslib-1.15/;/root/include/" -DCMAKE_PREFIX_PATH="`realpath thirdparty`" .. &&\
    make

RUN mkdir ../mydir
WORKDIR /mydir
ENV PATH /build/moni-align/build:$PATH