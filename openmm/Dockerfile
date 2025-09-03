# Base image
FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04

# Set working directory
WORKDIR /home

# Base packages
RUN apt-get update \
    && apt-get install -yq build-essential bzip2 wget csh \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Update conda
RUN conda update --quiet --yes conda

# OpenMM + tools
RUN . /opt/conda/etc/profile.d/conda.sh \
    && conda config --add channels conda-forge \
    && conda config --add channels omnia \
    && conda install -c conda-forge -c omnia \
        openmm cudatoolkit=12.2 pdbfixer parmed \
        mdtraj mdanalysis ambertools \
        numpy scipy 
    && conda clean -afy

CMD ["/bin/bash"]
