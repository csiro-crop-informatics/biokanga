FROM ubuntu:14.04
MAINTAINER Alex Whan <alex.whan@csiro.au>
RUN apt-get update
RUN apt-get install -y build-essential automake make
RUN mkdir /biokanga
ADD . /biokanga
WORKDIR "/biokanga"
RUN autoreconf -f -i
RUN ./configure
RUN make

