#Start with base Ubuntu image and install python
FROM ubuntu:24.04

#Add maintainer information
LABEL maintainer="Andrew Ramos" 
LABEL maintainer="ar39@sanger.ac.uk"

#Avoid promtps from apt 
ARG DEBIAN_FRONTEND=noninteractive

#Add image creation date 
ARG WHEN

# install basic dependencies, cleanup apt garbage.
RUN apt-get update \
    && apt-get upgrade -y \    
    && apt-get install -y \
        build-essential \
    && apt-get install -y python3-pip python3-full git-all \
    && apt-get -y install git-all \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

#Set working directory in the container
WORKDIR /usr/src/app

#Install SigProfilerAssignment package from PyPI
RUN pip3 install SigProfilerAssignment==0.2.3 --break-system-packages

#Test package installation and access 
RUN python3 -c 'import SigProfilerAssignment; print("SigProfilerAssignment package is ready!")'