## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin


RUN micromamba install -y -n base -c conda-forge -bioconda \
    conda-forge::r-tidyverse \
    conda-forge::datatable \
    conda-forge::r-optparse \
    conda-forge::r-janitor \
    bioconda::r-bedr \
    conda-forge::r-arrow \
    conda-forge::r-rfast \
    bioconda::r-optparse

ENV R_ENV=r-env 

RUN Rscript -e 'remotes::install_github("gabraham/plink2R@main")'

