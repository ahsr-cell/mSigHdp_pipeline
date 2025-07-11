FROM ubuntu:24.04

LABEL maintainer="Andrew Ramos" 
LABEL maintainer="ar39@sanger.ac.uk"

ARG DEBIAN_FRONTEND=noninteractive

ARG WHEN

# install basic dependencies, cleanup apt garbage.
RUN apt-get update \
    && apt-get upgrade -y \    
    && apt-get install -y \
        build-essential \
        curl \
        git \
        wget \
        default-jdk \
        libxml2 \
        libxml2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libssl-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libgit2-dev \
        r-base \
        r-cran-devtools \
        r-cran-biocmanager \
        r-cran-tidyverse \
        r-cran-reticulate \
        r-cran-rcurl \
        r-cran-rcpparmadillo \
        r-cran-xml \
        r-cran-xml2 \
        r-cran-git2r \
        r-cran-gert \
        r-cran-openssl \ 
        r-cran-textshaping \
        r-cran-systemfonts \
        r-cran-argparse \
        r-cran-r6 \
        r-cran-findpython \
        r-cran-jsonlite \
        r-bioc-rhtslib \
        r-bioc-genomeinfodb \
        r-bioc-biostrings \
        r-bioc-genomicranges \
        r-bioc-summarizedexperiment \
        r-bioc-rsamtools \
        r-bioc-genomicalignments \
        r-bioc-rtracklayer \
        r-bioc-bsgenome \
    && apt-get install -y python3-dev pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN R -e "devtools::install_git('https://github.com/ahsr-cell/hdpx'); \
    options(repos = c(CRAN='https://cloud.r-project.org')); \
    devtools::install_github('Rozen-Lab/mSigTools@v1.0.7'); \
    devtools::install_github('steverozen/mSigHdp@v2.1.2-bis', dependencies=FALSE)"

RUN R -e "library(dplyr); \
    library(mSigHdp); \ 
    print('All packages are ready to be used!')"
