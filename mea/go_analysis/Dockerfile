FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3 python3-pip libfontconfig1-dev libcurl4-openssl-dev libssl-dev gfortran libblas-dev liblapack-dev

RUN mkdir -p /opt/software/setup/R
ADD install_packages.R /opt/software/setup/R/
RUN Rscript /opt/software/setup/R/install_packages.R

WORKDIR /app

COPY ORA_redundancy_batch.R /app/

# add /app/ to path
