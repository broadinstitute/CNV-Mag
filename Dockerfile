# Use base image with miniconda3 installed
FROM continuumio/miniconda3
LABEL org.opencontainers.image.authors="gaoyueya@broadinstitute.org"

# Specify Workdir
WORKDIR /BaseImage

# Install bedtools
RUN conda install -y -c bioconda bioconda/label/main::bedtools && conda clean --all


# Create the environment
COPY Mag_env.yml .
RUN conda env create -f Mag_env.yml


# Copy the Cov Viz scripts to the container
RUN mkdir /BaseImage/CNV-Mag
COPY MagScripts/*py /BaseImage/CNV-Mag/
COPY MagRef/* /BaseImage/MagRef/