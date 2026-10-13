FROM ubuntu:24.04 AS builder

USER root

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends ca-certificates
RUN apt-get install -yq --no-install-recommends wget curl

RUN apt-get install -yq --no-install-recommends cmake

RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends pkg-config
RUN apt-get install -yq --no-install-recommends gcc g++

RUN apt-get install -yq --no-install-recommends r-base r-base-dev

RUN apt-get install -yq --no-install-recommends zlib1g-dev
RUN apt-get install -yq --no-install-recommends libbz2-dev
RUN apt-get install -yq --no-install-recommends liblzma-dev
RUN apt-get install -yq --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -yq --no-install-recommends libncurses-dev
RUN apt-get install -yq --no-install-recommends libssl-dev
RUN apt-get install -yq --no-install-recommends libblas-dev
RUN apt-get install -yq --no-install-recommends liblapack-dev
RUN apt-get install -yq --no-install-recommends gfortran
RUN apt-get install -yq --no-install-recommends libxml2-dev
RUN apt-get install -yq --no-install-recommends libgsl-dev
RUN apt-get install -yq --no-install-recommends libperl-dev
RUN apt-get install -yq --no-install-recommends libpng-dev

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT=/opt/wtsi-cgp
ENV PATH=$OPT/bin:$PATH
ENV R_LIBS=$OPT/R-lib
ENV R_LIBS_USER=$R_LIBS
ENV LD_LIBRARY_PATH=$OPT/lib
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# build tools from other repos
ADD versions.sh ./
ADD setup-fn.sh ./
ADD build-scripts/libInstall.R build-scripts/
ADD build-scripts/build-external-tools.sh build-scripts/
RUN bash build-scripts/build-external-tools.sh $OPT

# Install deepSNV
RUN mkdir -p "/opt/wtsi-cgp/R-lib"

RUN Rscript -e 'install.packages(c("remotes", "BiocManager"))'

RUN Rscript -e 'library("remotes"); remotes::install_version("curl", "5.2.1", lib = "/opt/wtsi-cgp/R-lib", lib.loc = "/opt/wtsi-cgp/R-lib")'
RUN Rscript -e 'library("remotes"); remotes::install_version("httr", "1.4.7", lib = "/opt/wtsi-cgp/R-lib", lib.loc = "/opt/wtsi-cgp/R-lib")'
RUN Rscript -e 'library("BiocManager"); BiocManager::install("VGAM", version = "3.18", update = FALSE, lib = "/opt/wtsi-cgp/R-lib",  lib.loc = "/opt/wtsi-cgp/R-lib")'
RUN Rscript -e 'library("BiocManager"); BiocManager::install("deepSNV", version = "3.18", update = FALSE, lib = "/opt/wtsi-cgp/R-lib",  lib.loc = "/opt/wtsi-cgp/R-lib")'
RUN Rscript -e 'library("BiocManager"); BiocManager::install("vcfR", version = "3.18", update = FALSE, lib = "/opt/wtsi-cgp/R-lib",  lib.loc = "/opt/wtsi-cgp/R-lib")'

# build the tools in this repo, separate to reduce build time on errors
COPY . .
ADD build-scripts/build-local.sh build-scripts/
RUN bash build-scripts/build-local.sh $OPT

FROM ubuntu:24.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="1.0.1" \
      description="nanoseq docker"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
curl \
wget \
make \
g++ \
gcc \
gfortran \
libblas-dev \
liblapack-dev \
ca-certificates \
time \
zlib1g \
libz-dev \
python3 \
libxml2 \
libgsl27 \
libperl5.38t64 \
libcapture-tiny-perl \
libfile-which-perl \
libpng16-16 \
parallel \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN apt-get install -yq --no-install-recommends r-base

ADD build-scripts/libInstall2.R build-scripts/
RUN Rscript build-scripts/libInstall2.R

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT=/opt/wtsi-cgp
ENV PATH=$OPT/bin:$PATH
ENV R_LIBS=$OPT/R-lib
ENV R_LIBS_USER=$R_LIBS
ENV LD_LIBRARY_PATH=$OPT/lib
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

CMD ["/bin/bash"]
