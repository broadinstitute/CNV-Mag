# Use base image with miniconda3 installed
FROM continuumio/miniconda3
LABEL org.opencontainers.image.authors="gaoyueya@broadinstitute.org"

# Specify Workdir
WORKDIR /BaseImage

# Install bedtools
RUN conda install -y -c bioconda bioconda/label/main::bedtools && conda clean --all

# Install gcloud CLI
# This a copy and paste from Galuoises' answer in stackoverflow
# Reference: https://stackoverflow.com/questions/28372328/how-to-install-the-google-cloud-sdk-in-a-docker-image
RUN apt-get update && \
    apt-get install -y curl gnupg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && \
    apt-get install google-cloud-sdk -y

# Create the environment
COPY Mag_env.yml .
RUN conda env create -f Mag_env.yml


# Copy the Cov Viz scripts to the container
RUN mkdir /BaseImage/CNV-Mag
COPY MagScripts/*py /BaseImage/CNV-Mag/
COPY MagRef/* /BaseImage/MagRef/