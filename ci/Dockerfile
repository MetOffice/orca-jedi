# Container image that runs your code
FROM jcsda/docker-gnu-openmpi-dev:revert AS base

RUN apt-get update \
    && apt-get -y install libhdf5-serial-dev \
    && apt-get -y install libnetcdf-dev libnetcdff-dev libnetcdf-c++4-dev \
    && cd /root \
    && git clone https://github.com/lz4/lz4.git \
    && cd lz4 \
    && git checkout tags/v1.9.3 \
    && make \
    && make --eval prefix=/opt install \
    && export LZ4_ROOT=/opt
