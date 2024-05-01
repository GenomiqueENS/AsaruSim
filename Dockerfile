# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

LABEL maintainer="Ali Hamraoui <hamraoui@bio.ens.psl.eu>"

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install --yes \
        build-essential \
        libncursesw5-dev \
        libgdbm-dev \
        libnss3-dev \
        libssl-dev \
        libsqlite3-dev \
        libreadline-dev \
        libffi-dev \
        libbz2-dev \
        wget \
        curl \
        default-jre \
        git \
        r-base && \
    wget https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz && \
    tar -xzf Python-3.11.0.tgz && \
    cd Python-3.11.0 && \
    ./configure --enable-optimizations && \
    make -j 8 && \
    make altinstall && \
    cd .. && \
    rm -rf Python-3.11.0.tgz Python-3.11.0


RUN python3.11 -m pip install --upgrade pip && \
    pip3 install numpy pandas tqdm biopython pyfaidx scipy plotly && \
    git clone https://github.com/rrwick/Badread.git /tmp/Badread && \
    python3.11 -m pip install /tmp/Badread && \
    rm -rf /tmp/Badread && \
    curl -fsSL get.nextflow.io | bash && \
    mv nextflow /usr/local/bin && \
    apt remove --yes git && \
    apt autoremove --yes && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    Rscript -e "install.packages(c('devtools', 'remotes', 'dplyr'), repos='https://cloud.r-project.org')" && \
    Rscript -e "remotes::install_gitlab('sysbiobig/sparsim', build_opts = c('--no-resave-data', '--no-manual'), build_vignettes = FALSE)"
