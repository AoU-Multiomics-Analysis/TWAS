## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin


RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::datatable \
    conda-forge::r-bedr \
    conda-forge::r-arrow \
    conda-forge::r-optparse \
    conda-forge::r-rfast=2.1.5.2

ENV R_ENV=r-env 
COPY scripts/*  . 
#RUN Rscript -e 'remotes::install_github("gabraham/plink2R@main")'
#RUN Rscript -e 'install.packages("Rfast)'

