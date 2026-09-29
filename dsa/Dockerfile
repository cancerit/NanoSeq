# Builder stage
FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    software-properties-common \
    ca-certificates \
    lsb-release \
    && rm -rf /var/lib/apt/lists/*
RUN wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-20 main" > /etc/apt/sources.list.d/llvm.list
RUN apt-get update && apt-get install -y \
    clang-20 \
    lldb-20 \
    lld-20 \
    libc++-20-dev \
    libc++abi-20-dev \
    libclang-20-dev \
    libclang-common-20-dev \
    && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang /usr/bin/clang-20 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-20 100 \
    && update-alternatives --install /usr/bin/lld lld /usr/bin/lld-20 100
RUN apt-get update && apt-get install -y \
    libhts-dev \
    libdeflate-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y \
    make \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
COPY CMakeLists.txt ./
COPY src/ ./src
RUN mkdir build \
    && cd build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release \
    && cmake --build .

# Runtime stage
FROM ubuntu:24.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    ca-certificates \
    && wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-20 main" > /etc/apt/sources.list.d/llvm.list \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y --no-install-recommends \
    libc++1-20 \
    libc++abi1-20 \
    libhts3 \
    libdeflate0 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libssl3 \
    strace \
    && apt-get purge -y wget gnupg \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
COPY --from=builder /workspace/build/dsa ./dsa
USER ubuntu
ENTRYPOINT ["/workspace/dsa"]

# Runtime stage with Python wrapper
FROM runtime AS runtime-python
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends python3.14 \
    && apt-get purge -y software-properties-common \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3.14 /usr/bin/python3

COPY scripts/ /workspace/scripts/

USER ubuntu

# Expose scripts
ENV PATH="/workspace/scripts:$PATH"

# Required by run_dsa.py
ENV DSA_EXE="/workspace/dsa"

CMD ["run_dsa.py"]
