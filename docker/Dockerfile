FROM ubuntu:20.04

RUN apt-get update -y && \
  apt-get install -y gcc libblas-dev liblapacke-dev make vim

WORKDIR /app
COPY . .

ENV HOSTNAME=docker
RUN mkdir bin
# RUN make all && \
#   make run