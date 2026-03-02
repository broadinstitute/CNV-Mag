# Use base image with micromamba installed
FROM mambaorg/micromamba:2.5-amazon2023
LABEL org.opencontainers.image.authors="gaoyueya@broadinstitute.org"

# Install bedtools and bcftools
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    bedtools bcftools \
    && micromamba clean --all --yes

# Create the environment
COPY Mag_env.yml .
RUN micromamba create -n CNV-Mag -f Mag_env.yml

# Specify Workdir
WORKDIR /BaseImage

# Copy the Cov Viz scripts to the container
RUN mkdir /BaseImage/CNV-Mag
COPY MagScripts/*py /BaseImage/CNV-Mag/
COPY MagRef/* /BaseImage/MagRef/