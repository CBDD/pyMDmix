FROM ggutierrezbio/ambertools20
ARG CONDA_PREFIX="/opt/SOFT/miniconda3"
ENV PATH="$CONDA_PREFIX/bin:${PATH}"
ARG PATH="$CONDA_PREFIX/bin:${PATH}"

SHELL [ "/bin/bash", "--login", "-c" ]
RUN apt-get update && apt-get install -y wget git && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_PREFIX \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

COPY ./environment_p27.yml ./environment_p27.yml
RUN conda env create -f environment_p27.yml
RUN echo ". $CONDA_PREFIX/etc/profile.d/conda.sh" >> ~/.bash_profile
RUN echo "conda activate mdmix-env" >> ~/.bash_profile
WORKDIR /mnt
SHELL ["conda", "run", "-n", "mdmix-env", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "-n", "mdmix-env", "mdmix"]
