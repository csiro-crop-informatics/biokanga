FROM ubuntu:16.04
MAINTAINER Rad Suchecki <rad.suchecki@csiro.au>
RUN apt-get update && apt-get install -y wget build-essential automake make

#overwrite at build time e.g. with --build-arg version=4.3.11
ARG VERSION=4.3.10

#Install biokanga
WORKDIR /biokanga
RUN wget https://github.com/csiro-crop-informatics/biokanga/archive/v${VERSION}.tar.gz \
  && tar xzf v${VERSION}.tar.gz \
  && cd biokanga-${VERSION} \
  && autoreconf -f -i \
  && ./configure \
  && make
ENV PATH "$PATH:/biokanga/biokanga-${VERSION}/biokanga"