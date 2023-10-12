FROM  ubuntu:18.04 as builder

USER  root

# ALL tool versions used by opt-build.sh
ENV VER_SAMTOOLS="1.18"
ENV VER_HTSLIB="1.18"
ENV VER_BCFTOOLS="1.18"
ENV VER_VERIFYBAMID="2.0.1"
ENV VER_LIBDEFLATE="v1.18"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends g++
RUN apt-get install -yq --no-install-recommends ca-certificates
RUN apt-get install -yq --no-install-recommends wget

# install latest cmake so opt-build.sh works - the initial installs will also help install R
RUN apt-get install -yq --no-install-recommends software-properties-common lsb-release
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
RUN apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main"
RUN apt-get install -yq --no-install-recommends cmake=3.25.2-0kitware1ubuntu18.04.1

RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends pkg-config

# if ubuntu 18.04
RUN apt install -yq --no-install-recommends dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get install -yq --no-install-recommends r-recommended=4.1.3-1.1804.0
RUN apt-get install -yq --no-install-recommends r-base=4.1.3-1.1804.0
# if ubuntu 22.04
# RUN apt-get install -yq --no-install-recommends r-base=4.1.2-1ubuntu2

RUN apt-get install -yq --no-install-recommends zlib1g-dev
RUN apt-get install -yq --no-install-recommends libbz2-dev
RUN apt-get install -yq --no-install-recommends liblzma-dev
RUN apt-get install -yq --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -yq --no-install-recommends libncurses5-dev
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

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# build tools from other repos
ADD build/libInstall.R build/
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:18.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="1.0.1" \
      description="nanoseq docker"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
wget \
ca-certificates \
time \
zlib1g \
python3 \
libxml2 \
libgsl23 \
libperl5.26 \
libcapture-tiny-perl \
libfile-which-perl \
libpng16-16 \
parallel \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN apt install -yq software-properties-common --no-install-recommends dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get install -yq --no-install-recommends r-recommended=4.1.3-1.1804.0
RUN apt-get install -yq --no-install-recommends r-base=4.1.3-1.1804.0
RUN apt-get install -yq --no-install-recommends r-cran-ggplot2
RUN apt-get install -yq --no-install-recommends r-cran-data.table
RUN apt-get install -yq --no-install-recommends r-cran-epitools
RUN apt-get install -yq --no-install-recommends r-cran-gridextra
RUN apt-get install -yq --no-install-recommends r-cran-seqinr

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
