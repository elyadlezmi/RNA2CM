# Set the base image to debian
FROM ubuntu:focal

# File Author / Maintainer
MAINTAINER Elyad Lezmi <elyad.lezmi@mail.huji.ac.il>

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    unzip \
    wget \
    samtools \
    tabix \
    r-base \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-openssl-dev \ 
    libssl-dev \
    libxml2-dev \
    rna-star \
    bcftools \
    python \
    python3-pip \
    openjdk-8-jdk

RUN pip3 install pandas

RUN R -q -e "install.packages('BiocManager', repos='https://cran.r-project.org')"
RUN R -q -e 'BiocManager::install(c("Rsamtools", "GenomicAlignments", "BiocParallel", "futile.logger"))'
RUN wget https://github.com/PeeperLab/XenofilteR/releases/download/v1.6/XenofilteR_1.6.tar.gz
RUN R CMD INSTALL XenofilteR_1.6.tar.gz
RUN rm XenofilteR_1.6.tar.gz

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
RUN unzip Trimmomatic-0.39.zip
RUN rm Trimmomatic-0.39.zip

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
RUN unzip gatk-4.1.3.0.zip
RUN rm gatk-4.1.3.0.zip
