FROM ubuntu:16.04

RUN apt-get update && apt-get install -y build-essential automake make

RUN mkdir /biokanga
ADD . /biokanga
WORKDIR /biokanga
RUN autoreconf -f -i \
  && ./configure \
  && make

RUN apt remove --purge --yes \
    make \
    automake \
    build-essential \
  && apt autoremove --purge --yes

ENV PATH "/biokanga/biokanga:${PATH}"
