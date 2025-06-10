FROM continuumio/miniconda3 AS base

ARG PYTHON_VERSION=3.12

# Create environment and install ambertools
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# RUN apt update && apt install libnetcdf-dev -y

RUN conda install -y python=${PYTHON_VERSION} ambertools

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

CMD ["bash"]

FROM base AS development

COPY requirements-dev.txt .
RUN python -m pip install -r requirements-dev.txt
